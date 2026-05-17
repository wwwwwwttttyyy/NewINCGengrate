#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;

static constexpr double PI = 3.141592653589793238462643383279502884;

struct Corner {
    int v = 0;
    int sx = 0;
    int sy = 0;
};

struct Face {
    Corner c[3];
};

struct EdgeAdj {
    int a = 0;
    int b = 0;
    int shift_x = 0;
    int shift_y = 0;
    int f0 = -1;
    int e0 = -1;
    int f1 = -1;
    int e1 = -1;
};

struct Mesh {
    int nx = 0;
    int ny = 0;
    int N = 0;
    vector<Face> faces;
    vector<EdgeAdj> edges;
    vector<array<int, 3>> face_edges;
    vector<vector<int>> incident_faces;
};

struct State {
    vector<double> u;
    vector<double> r;
    vector<double> K;
    double E_K = 0.0;
    double max_abs_K = 0.0;
    double r_min = 0.0;
    double r_max = 0.0;
};

struct StatRow {
    int iter = 0;
    double E_K = 0.0;
    double max_abs_K = 0.0;
    double step = 0.0;
    double r_min = 0.0;
    double r_max = 0.0;
    int accepted = 1;
};

struct RicciResult {
    State state;
    vector<StatRow> stats;
    int iterations = 0;
    string stopped_reason;
    bool monotonic_E = true;
    int boundary_rejections = 0;
    double initial_E_K = 0.0;
    double initial_max_abs_K = 0.0;
    double u_min_actual = 0.0;
    double u_max_actual = 0.0;
    double radius_ratio = 0.0;
};

struct DevelopedFace {
    bool done = false;
    double x[3] = {};
    double y[3] = {};
    int sx[3] = {};
    int sy[3] = {};
};

struct DevelopResult {
    bool ran = false;
    int faces_count = 0;
    string status = "not requested";
    vector<DevelopedFace> faces;
    vector<double> theta_vertex;
    vector<int> theta_vertex_count;
    double global_vertex_max_spread = numeric_limits<double>::quiet_NaN();
    double global_vertex_rms_spread = numeric_limits<double>::quiet_NaN();
    double period_fit_rms = numeric_limits<double>::quiet_NaN();
    double period_A_x = numeric_limits<double>::quiet_NaN();
    double period_A_y = numeric_limits<double>::quiet_NaN();
    double period_B_x = numeric_limits<double>::quiet_NaN();
    double period_B_y = numeric_limits<double>::quiet_NaN();
    double theta_mean_abs = numeric_limits<double>::quiet_NaN();
    double theta_rms = numeric_limits<double>::quiet_NaN();
    double theta_max_abs = numeric_limits<double>::quiet_NaN();
};

struct Timings {
    double topology = 0.0;
    double ricci = 0.0;
    double develop = 0.0;
    double embed = 0.0;
    double total = 0.0;
};

struct Config {
    int nx = 16;
    int ny = 16;
    int flips = 0;
    unsigned long long seed = 1;
    int max_iter = 200000;
    double tol = 1e-10;
    string method = "ricci";
    string out = "out";
    bool test = false;

    bool bounded = false;
    double u_min = -2.0;
    double u_max = 2.0;

    bool develop = false;
    bool embed = false;

    int progress_every = 1000;
    int validate_every = 0;
    int stagnation_window = 5000;
    double stagnation_rel = 1e-12;

    bool scan_flips = false;
    vector<int> scan_list;
    vector<unsigned long long> scan_seeds;
    string scan_out = "scan.csv";
    bool scan_write_details = false;
    int threads = 4;

    int write_develop_faces = -1;
    int write_theta_vertex = -1;
    int write_coords = -1;
};

struct RunOutput {
    Config cfg;
    Mesh mesh;
    RicciResult ricci;
    DevelopResult develop;
    Timings time;
    int flips_accepted = 0;
    bool valid = true;
    string error;
};

struct ScanRow {
    int nx = 0, ny = 0, N = 0, E = 0, F = 0;
    int flips_requested = 0;
    int flips_accepted = 0;
    unsigned long long seed = 0;
    int ricci_iters = 0;
    double final_E_K = numeric_limits<double>::quiet_NaN();
    double final_max_abs_K = numeric_limits<double>::quiet_NaN();
    double radius_ratio = numeric_limits<double>::quiet_NaN();
    double develop_global_vertex_max_spread = numeric_limits<double>::quiet_NaN();
    double period_fit_rms = numeric_limits<double>::quiet_NaN();
    double develop_theta_mean_abs = numeric_limits<double>::quiet_NaN();
    double time_total_sec = 0.0;
    double time_ricci_sec = 0.0;
    double time_develop_sec = 0.0;
    string stopped_reason;
    string status = "ok";
};

static double now_sec() {
    using clock = chrono::steady_clock;
    static const auto t0 = clock::now();
    return chrono::duration<double>(clock::now() - t0).count();
}

static int imod(int x, int n) {
    int r = x % n;
    return r < 0 ? r + n : r;
}

static int vertex_id(int i, int j, int nx, int ny) {
    return imod(i, nx) + nx * imod(j, ny);
}

static int vertex_ix(int v, int nx) {
    return v % nx;
}

static int vertex_iy(int v, int nx) {
    return v / nx;
}

static int floor_div_periodic(int x, int n) {
    int base = imod(x, n);
    return (x - base) / n;
}

static Corner corner_from_unwrapped(int i, int j, int nx, int ny) {
    return {vertex_id(i, j, nx, ny), floor_div_periodic(i, nx), floor_div_periodic(j, ny)};
}

static Face make_face(Corner a, Corner b, Corner c) {
    Face f;
    f.c[0] = a;
    f.c[1] = b;
    f.c[2] = c;
    return f;
}

static uint64_t pack_edge_key(int a, int b) {
    if (a > b) swap(a, b);
    return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
}

static double clamp_unit(double x) {
    return max(-1.0, min(1.0, x));
}

static Mesh generate_initial_mesh(int nx, int ny) {
    Mesh mesh;
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.N = nx * ny;
    mesh.faces.reserve(size_t(2 * mesh.N));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Corner A = corner_from_unwrapped(i, j, nx, ny);
            Corner B = corner_from_unwrapped(i + 1, j, nx, ny);
            Corner C = corner_from_unwrapped(i, j + 1, nx, ny);
            Corner D = corner_from_unwrapped(i + 1, j + 1, nx, ny);
            mesh.faces.push_back(make_face(A, B, C));
            mesh.faces.push_back(make_face(B, D, C));
        }
    }
    return mesh;
}

static bool canonical_shift_for_edge(const Face& f, int e, int* sx, int* sy) {
    int p = e;
    int q = (e + 1) % 3;
    const Corner& a = f.c[p];
    const Corner& b = f.c[q];
    int lo = min(a.v, b.v);
    int hi = max(a.v, b.v);
    int dx = b.sx - a.sx;
    int dy = b.sy - a.sy;
    if (a.v == lo && b.v == hi) {
        *sx = dx;
        *sy = dy;
    } else if (a.v == hi && b.v == lo) {
        *sx = -dx;
        *sy = -dy;
    } else {
        return false;
    }
    return true;
}

static bool build_adjacency(Mesh& mesh, string* error = nullptr,
                            unordered_map<uint64_t, int>* edge_index_out = nullptr) {
    const int F = static_cast<int>(mesh.faces.size());
    mesh.edges.clear();
    mesh.face_edges.assign(F, {-1, -1, -1});
    mesh.incident_faces.assign(mesh.N, {});

    unordered_map<uint64_t, int> edge_index;
    edge_index.reserve(size_t(3 * F * 2 + 1));
    mesh.edges.reserve(size_t(3 * mesh.N));

    auto fail = [&](const string& msg) {
        if (error) *error = msg;
        return false;
    };

    for (int fi = 0; fi < F; ++fi) {
        const Face& f = mesh.faces[fi];
        for (int q = 0; q < 3; ++q) {
            if (f.c[q].v < 0 || f.c[q].v >= mesh.N) return fail("face has out-of-range vertex");
        }
        if (f.c[0].v == f.c[1].v || f.c[1].v == f.c[2].v || f.c[2].v == f.c[0].v) {
            return fail("face has duplicate vertex ids");
        }
        for (int q = 0; q < 3; ++q) mesh.incident_faces[f.c[q].v].push_back(fi);

        for (int e = 0; e < 3; ++e) {
            int a = f.c[e].v;
            int b = f.c[(e + 1) % 3].v;
            uint64_t key = pack_edge_key(a, b);
            int sx = 0, sy = 0;
            if (!canonical_shift_for_edge(f, e, &sx, &sy)) return fail("bad canonical shift");

            auto it = edge_index.find(key);
            if (it == edge_index.end()) {
                EdgeAdj adj;
                adj.a = min(a, b);
                adj.b = max(a, b);
                adj.shift_x = sx;
                adj.shift_y = sy;
                adj.f0 = fi;
                adj.e0 = e;
                int idx = static_cast<int>(mesh.edges.size());
                mesh.edges.push_back(adj);
                edge_index.emplace(key, idx);
                mesh.face_edges[fi][e] = idx;
            } else {
                EdgeAdj& adj = mesh.edges[it->second];
                if (adj.shift_x != sx || adj.shift_y != sy) return fail("inconsistent periodic edge shift");
                if (adj.f1 != -1) return fail("edge has more than two incident faces");
                adj.f1 = fi;
                adj.e1 = e;
                mesh.face_edges[fi][e] = it->second;
            }
        }
    }

    for (const EdgeAdj& e : mesh.edges) {
        if (e.f0 < 0 || e.f1 < 0) return fail("edge does not have exactly two incident faces");
    }

    if (edge_index_out) *edge_index_out = std::move(edge_index);
    return true;
}

static bool validate_mesh(Mesh& mesh, string* error = nullptr) {
    if (!build_adjacency(mesh, error)) return false;
    int V = mesh.N;
    int E = static_cast<int>(mesh.edges.size());
    int F = static_cast<int>(mesh.faces.size());
    if (V - E + F != 0) {
        if (error) *error = "Euler characteristic is not zero";
        return false;
    }
    if (F != 2 * V) {
        if (error) *error = "F != 2V";
        return false;
    }
    if (E != 3 * V) {
        if (error) *error = "E != 3V";
        return false;
    }
    return true;
}

static int corner_index_for_vertex(const Face& f, int v) {
    for (int i = 0; i < 3; ++i) {
        if (f.c[i].v == v) return i;
    }
    return -1;
}

static int opposite_index(const Face& f, int a, int b) {
    for (int i = 0; i < 3; ++i) {
        if (f.c[i].v != a && f.c[i].v != b) return i;
    }
    return -1;
}

static Corner shifted_corner(Corner c, int ox, int oy) {
    c.sx += ox;
    c.sy += oy;
    return c;
}

static bool attempt_random_flip(Mesh& mesh, mt19937_64& rng) {
    unordered_map<uint64_t, int> edge_index;
    string error;
    if (!build_adjacency(mesh, &error, &edge_index)) return false;
    if (mesh.edges.empty()) return false;

    uniform_int_distribution<int> pick(0, static_cast<int>(mesh.edges.size()) - 1);
    const EdgeAdj edge = mesh.edges[pick(rng)];
    if (edge.f0 < 0 || edge.f1 < 0) return false;

    const Face& face1 = mesh.faces[edge.f0];
    const Face& face2 = mesh.faces[edge.f1];
    int ia1 = corner_index_for_vertex(face1, edge.a);
    int ib1 = corner_index_for_vertex(face1, edge.b);
    int ia2 = corner_index_for_vertex(face2, edge.a);
    int ib2 = corner_index_for_vertex(face2, edge.b);
    int io1 = opposite_index(face1, edge.a, edge.b);
    int io2 = opposite_index(face2, edge.a, edge.b);
    if (ia1 < 0 || ib1 < 0 || ia2 < 0 || ib2 < 0 || io1 < 0 || io2 < 0) return false;

    Corner a = face1.c[ia1];
    Corner b = face1.c[ib1];
    Corner c = face1.c[io1];
    int ox = a.sx - face2.c[ia2].sx;
    int oy = a.sy - face2.c[ia2].sy;
    Corner b2 = shifted_corner(face2.c[ib2], ox, oy);
    if (b2.sx != b.sx || b2.sy != b.sy) return false;
    Corner d = shifted_corner(face2.c[io2], ox, oy);
    if (c.v == d.v) return false;
    if (edge_index.find(pack_edge_key(c.v, d.v)) != edge_index.end()) return false;

    mesh.faces[edge.f0] = make_face(c, d, a);
    mesh.faces[edge.f1] = make_face(d, c, b);
    return true;
}

static int randomize_topology(Mesh& mesh, int flips, int validate_every, mt19937_64& rng) {
    int accepted = 0;
    for (int i = 0; i < flips; ++i) {
        if (attempt_random_flip(mesh, rng)) {
            accepted++;
            if (validate_every > 0 && accepted % validate_every == 0) {
                string error;
                if (!validate_mesh(mesh, &error)) throw runtime_error("validation failed during flips: " + error);
            }
        }
    }
    return accepted;
}

static void normalize_u(vector<double>& u) {
    double mean = accumulate(u.begin(), u.end(), 0.0) / max<size_t>(u.size(), 1);
    for (double& x : u) x -= mean;
}

static vector<double> random_log_radii(int N, mt19937_64& rng, double stddev = 0.1) {
    normal_distribution<double> normal(0.0, stddev);
    vector<double> u(N);
    for (double& x : u) x = normal(rng);
    normalize_u(u);
    return u;
}

static double safe_exp(double x) {
    if (x > 700.0) x = 700.0;
    if (x < -700.0) x = -700.0;
    return exp(x);
}

static array<double, 3> angles_from_radii(double ra, double rb, double rc) {
    auto angle_at = [](double ri, double rj, double rk) {
        double a = ri + rj;
        double b = ri + rk;
        double c = rj + rk;
        double scale = max(a, max(b, c));
        if (scale <= 0.0 || !isfinite(scale)) return numeric_limits<double>::quiet_NaN();
        a /= scale;
        b /= scale;
        c /= scale;
        double denom = 2.0 * a * b;
        if (denom <= 0.0) return numeric_limits<double>::quiet_NaN();
        return acos(clamp_unit((a * a + b * b - c * c) / denom));
    };
    return {angle_at(ra, rb, rc), angle_at(rb, rc, ra), angle_at(rc, ra, rb)};
}

static void compute_state_inplace(const Mesh& mesh, const vector<double>& u, State& s) {
    const int N = mesh.N;
    s.u = u;
    s.r.assign(N, 0.0);
    s.K.assign(N, 2.0 * PI);
    s.r_min = numeric_limits<double>::infinity();
    s.r_max = 0.0;

    for (int i = 0; i < N; ++i) {
        s.r[i] = safe_exp(u[i]);
        s.r_min = min(s.r_min, s.r[i]);
        s.r_max = max(s.r_max, s.r[i]);
    }

    for (const Face& f : mesh.faces) {
        int a = f.c[0].v;
        int b = f.c[1].v;
        int c = f.c[2].v;
        auto ang = angles_from_radii(s.r[a], s.r[b], s.r[c]);
        s.K[a] -= ang[0];
        s.K[b] -= ang[1];
        s.K[c] -= ang[2];
    }

    s.E_K = 0.0;
    s.max_abs_K = 0.0;
    for (double k : s.K) {
        s.E_K += 0.5 * k * k;
        s.max_abs_K = max(s.max_abs_K, fabs(k));
    }
}

static bool inside_bounds(const vector<double>& u, double lo, double hi) {
    for (double x : u) {
        if (x < lo || x > hi) return false;
    }
    return true;
}

static StatRow make_stat(int iter, const State& s, double step, int accepted) {
    return {iter, s.E_K, s.max_abs_K, step, s.r_min, s.r_max, accepted};
}

static void finalize_ricci(RicciResult& rr) {
    if (!rr.state.u.empty()) {
        auto mm = minmax_element(rr.state.u.begin(), rr.state.u.end());
        rr.u_min_actual = *mm.first;
        rr.u_max_actual = *mm.second;
    }
    rr.radius_ratio = rr.state.r_min > 0.0 ? rr.state.r_max / rr.state.r_min
                                            : numeric_limits<double>::infinity();
}

static RicciResult run_ricci(const Mesh& mesh, vector<double> u, const Config& cfg, bool keep_stats, bool verbose) {
    normalize_u(u);
    RicciResult rr;
    State current, candidate;
    compute_state_inplace(mesh, u, current);
    rr.initial_E_K = current.E_K;
    rr.initial_max_abs_K = current.max_abs_K;
    if (keep_stats) rr.stats.push_back(make_stat(0, current, 0.01, 1));

    double dt = 0.01;
    const double min_dt = 1e-14;
    vector<double> cand_u(mesh.N, 0.0);
    double window_E = current.E_K;
    int window_start = 0;

    if (verbose && cfg.progress_every > 0) {
        cout << scientific << setprecision(6)
             << "iter 0 E_K " << current.E_K
             << " max_abs_K " << current.max_abs_K
             << " step " << dt
             << " r_min " << current.r_min
             << " r_max " << current.r_max << '\n';
    }

    if (current.max_abs_K < cfg.tol) {
        rr.state = std::move(current);
        rr.stopped_reason = "converged: max_abs_K < tol at iter 0";
        finalize_ricci(rr);
        return rr;
    }

    for (int iter = 1; iter <= cfg.max_iter; ++iter) {
        bool accepted = false;
        double accepted_step = dt;

        while (!accepted) {
            if (dt < min_dt) {
                rr.state = std::move(current);
                rr.iterations = iter - 1;
                rr.stopped_reason = "stalled: dt < 1e-14";
                if (keep_stats) rr.stats.push_back(make_stat(rr.iterations, rr.state, dt, 0));
                finalize_ricci(rr);
                return rr;
            }

            for (int i = 0; i < mesh.N; ++i) cand_u[i] = u[i] - dt * current.K[i];
            normalize_u(cand_u);
            if (cfg.bounded && !inside_bounds(cand_u, cfg.u_min, cfg.u_max)) {
                rr.boundary_rejections++;
                dt *= 0.5;
                continue;
            }

            compute_state_inplace(mesh, cand_u, candidate);
            if (!isfinite(candidate.E_K) || !isfinite(candidate.max_abs_K)) {
                rr.state = std::move(current);
                rr.iterations = iter - 1;
                rr.stopped_reason = "failed: NaN or Inf in curvature";
                finalize_ricci(rr);
                return rr;
            }

            if (candidate.E_K < current.E_K) {
                accepted = true;
                accepted_step = dt;
                rr.monotonic_E = rr.monotonic_E && (candidate.E_K <= current.E_K + 1e-15);
                u.swap(cand_u);
                current = std::move(candidate);
                dt = min(dt * 1.05, 1.0);
            } else {
                dt *= 0.5;
            }
        }

        rr.iterations = iter;
        bool is_progress = cfg.progress_every > 0 && (iter % cfg.progress_every == 0);
        if (keep_stats && (is_progress || current.max_abs_K < cfg.tol)) {
            rr.stats.push_back(make_stat(iter, current, accepted_step, 1));
        }
        if (verbose && is_progress) {
            cout << scientific << setprecision(6)
                 << "iter " << iter
                 << " E_K " << current.E_K
                 << " max_abs_K " << current.max_abs_K
                 << " step " << accepted_step
                 << " r_min " << current.r_min
                 << " r_max " << current.r_max << '\n';
        }

        if (current.max_abs_K < cfg.tol) {
            rr.state = std::move(current);
            rr.stopped_reason = "converged: max_abs_K < tol";
            finalize_ricci(rr);
            return rr;
        }

        if (cfg.stagnation_window > 0 && iter - window_start >= cfg.stagnation_window) {
            double rel = (window_E - current.E_K) / max(window_E, 1e-300);
            if (rel < cfg.stagnation_rel) {
                rr.state = std::move(current);
                rr.stopped_reason = "stalled_by_stagnation";
                finalize_ricci(rr);
                return rr;
            }
            window_E = current.E_K;
            window_start = iter;
        }
    }

    rr.state = std::move(current);
    rr.stopped_reason = "reached max_iter";
    if (keep_stats && (rr.stats.empty() || rr.stats.back().iter != rr.iterations)) {
        rr.stats.push_back(make_stat(rr.iterations, rr.state, dt, 1));
    }
    finalize_ricci(rr);
    return rr;
}

static double target_length(const vector<double>& r, int a, int b) {
    return r[a] + r[b];
}

static double dist2(double ax, double ay, double bx, double by) {
    double dx = ax - bx;
    double dy = ay - by;
    return dx * dx + dy * dy;
}

static bool triangle_third_point(double ax, double ay, double bx, double by, double len_ac, double len_bc,
                                 double refx, double refy, double& cx, double& cy) {
    double ex = bx - ax;
    double ey = by - ay;
    double L = sqrt(ex * ex + ey * ey);
    if (L <= 1e-14) return false;
    double t = (len_ac * len_ac + L * L - len_bc * len_bc) / (2.0 * L);
    double h2 = len_ac * len_ac - t * t;
    if (h2 < 0.0 && h2 > -1e-10) h2 = 0.0;
    if (h2 < 0.0) return false;
    double h = sqrt(h2);
    double ux = ex / L;
    double uy = ey / L;
    double px = -uy;
    double py = ux;
    double c1x = ax + t * ux + h * px;
    double c1y = ay + t * uy + h * py;
    double c2x = ax + t * ux - h * px;
    double c2y = ay + t * uy - h * py;
    double side_ref = ex * (refy - ay) - ey * (refx - ax);
    double side1 = ex * (c1y - ay) - ey * (c1x - ax);
    if (fabs(side_ref) < 1e-14) {
        cx = c1x;
        cy = c1y;
    } else if (side_ref * side1 < 0.0) {
        cx = c1x;
        cy = c1y;
    } else {
        cx = c2x;
        cy = c2y;
    }
    return true;
}

struct LiftKey {
    int v = 0;
    int sx = 0;
    int sy = 0;
    bool operator==(const LiftKey& other) const {
        return v == other.v && sx == other.sx && sy == other.sy;
    }
};

struct LiftKeyHash {
    size_t operator()(const LiftKey& k) const {
        uint64_t h = uint64_t(uint32_t(k.v));
        h ^= uint64_t(uint32_t(k.sx + 1048576)) << 21;
        h ^= uint64_t(uint32_t(k.sy + 1048576)) << 42;
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        return size_t(h);
    }
};

struct SpreadAcc {
    int count = 0;
    double refx = 0.0;
    double refy = 0.0;
    double max_spread = 0.0;
    double sumsq = 0.0;
};

static double angle_between(double ax, double ay, double bx, double by) {
    double na = sqrt(ax * ax + ay * ay);
    double nb = sqrt(bx * bx + by * by);
    if (na <= 0.0 || nb <= 0.0) return numeric_limits<double>::quiet_NaN();
    return acos(clamp_unit((ax * bx + ay * by) / (na * nb)));
}

static DevelopResult develop_mesh(const Mesh& mesh, const vector<double>& r) {
    DevelopResult dr;
    dr.ran = true;
    dr.status = "ok";
    const int F = static_cast<int>(mesh.faces.size());
    dr.faces.assign(F, {});
    dr.theta_vertex.assign(mesh.N, 0.0);
    dr.theta_vertex_count.assign(mesh.N, 0);
    if (F == 0) {
        dr.status = "empty mesh";
        return dr;
    }

    queue<int> q;
    auto seed_face = [&](int fi) {
        const Face& f = mesh.faces[fi];
        int a = f.c[0].v, b = f.c[1].v, c = f.c[2].v;
        double lab = target_length(r, a, b);
        double lac = target_length(r, a, c);
        double lbc = target_length(r, b, c);
        double x2 = (lac * lac + lab * lab - lbc * lbc) / (2.0 * max(lab, 1e-14));
        double y2_sq = lac * lac - x2 * x2;
        if (y2_sq < 0.0 && y2_sq > -1e-10) y2_sq = 0.0;
        if (lab <= 0.0 || y2_sq < 0.0) return false;
        DevelopedFace& df = dr.faces[fi];
        df.done = true;
        df.x[0] = 0.0; df.y[0] = 0.0;
        df.x[1] = lab; df.y[1] = 0.0;
        df.x[2] = x2; df.y[2] = sqrt(y2_sq);
        for (int i = 0; i < 3; ++i) {
            df.sx[i] = f.c[i].sx;
            df.sy[i] = f.c[i].sy;
        }
        q.push(fi);
        dr.faces_count++;
        return true;
    };

    if (!seed_face(0)) {
        dr.status = "failed to seed first face";
        return dr;
    }

    while (!q.empty()) {
        int fi = q.front();
        q.pop();
        const Face& f = mesh.faces[fi];
        const DevelopedFace& df = dr.faces[fi];

        for (int e = 0; e < 3; ++e) {
            int edge_idx = mesh.face_edges[fi][e];
            if (edge_idx < 0) continue;
            const EdgeAdj& edge = mesh.edges[edge_idx];
            int nf = edge.f0 == fi ? edge.f1 : edge.f0;
            if (nf < 0 || dr.faces[nf].done) continue;

            int p = e;
            int rloc = (e + 1) % 3;
            int cur_opp = (e + 2) % 3;
            int vp = f.c[p].v;
            int vr = f.c[rloc].v;

            const Face& g = mesh.faces[nf];
            int np = corner_index_for_vertex(g, vp);
            int nr = corner_index_for_vertex(g, vr);
            int no = opposite_index(g, vp, vr);
            if (np < 0 || nr < 0 || no < 0) continue;

            int ox = df.sx[p] - g.c[np].sx;
            int oy = df.sy[p] - g.c[np].sy;
            Corner gr_shifted = shifted_corner(g.c[nr], ox, oy);
            if (gr_shifted.sx != df.sx[rloc] || gr_shifted.sy != df.sy[rloc]) continue;

            DevelopedFace nd;
            nd.done = true;
            nd.x[np] = df.x[p]; nd.y[np] = df.y[p];
            nd.x[nr] = df.x[rloc]; nd.y[nr] = df.y[rloc];
            nd.sx[np] = df.sx[p]; nd.sy[np] = df.sy[p];
            nd.sx[nr] = df.sx[rloc]; nd.sy[nr] = df.sy[rloc];

            Corner go = shifted_corner(g.c[no], ox, oy);
            nd.sx[no] = go.sx;
            nd.sy[no] = go.sy;
            double cx = 0.0, cy = 0.0;
            bool ok = triangle_third_point(df.x[p], df.y[p], df.x[rloc], df.y[rloc],
                                           target_length(r, vp, g.c[no].v),
                                           target_length(r, vr, g.c[no].v),
                                           df.x[cur_opp], df.y[cur_opp], cx, cy);
            if (!ok) continue;
            nd.x[no] = cx;
            nd.y[no] = cy;
            dr.faces[nf] = nd;
            dr.faces_count++;
            q.push(nf);
        }
    }

    if (dr.faces_count != F) {
        dr.status = "incomplete face development";
    }

    unordered_map<LiftKey, SpreadAcc, LiftKeyHash> spread;
    spread.reserve(size_t(3 * F * 2 + 1));
    double spread_sumsq = 0.0;
    int spread_count = 0;

    struct RefOcc {
        bool have = false;
        int sx = 0, sy = 0;
        double x = 0.0, y = 0.0;
    };
    vector<RefOcc> refs(mesh.N);
    double S00 = 0.0, S01 = 0.0, S11 = 0.0;
    double rx0 = 0.0, rx1 = 0.0, ry0 = 0.0, ry1 = 0.0;
    int constraints = 0;

    auto visit_occurrence = [&](int v, int sx, int sy, double x, double y) {
        LiftKey key{v, sx, sy};
        SpreadAcc& acc = spread[key];
        if (acc.count == 0) {
            acc.refx = x;
            acc.refy = y;
        } else {
            double d2 = dist2(x, y, acc.refx, acc.refy);
            double d = sqrt(d2);
            acc.max_spread = max(acc.max_spread, d);
            spread_sumsq += d2;
            spread_count++;
        }
        acc.count++;

        RefOcc& ref = refs[v];
        if (!ref.have) {
            ref.have = true;
            ref.sx = sx;
            ref.sy = sy;
            ref.x = x;
            ref.y = y;
        } else {
            double dsx = static_cast<double>(sx - ref.sx);
            double dsy = static_cast<double>(sy - ref.sy);
            if (dsx != 0.0 || dsy != 0.0) {
                double dx = x - ref.x;
                double dy = y - ref.y;
                S00 += dsx * dsx;
                S01 += dsx * dsy;
                S11 += dsy * dsy;
                rx0 += dsx * dx;
                rx1 += dsy * dx;
                ry0 += dsx * dy;
                ry1 += dsy * dy;
                constraints++;
            }
        }
    };

    double theta_sum = 0.0;
    double theta_sumsq = 0.0;
    int theta_count = 0;

    for (int fi = 0; fi < F; ++fi) {
        if (!dr.faces[fi].done) continue;
        const Face& f = mesh.faces[fi];
        const DevelopedFace& df = dr.faces[fi];
        for (int i = 0; i < 3; ++i) {
            visit_occurrence(f.c[i].v, df.sx[i], df.sy[i], df.x[i], df.y[i]);
        }

        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            double actual = angle_between(df.x[j] - df.x[i], df.y[j] - df.y[i],
                                          df.x[k] - df.x[i], df.y[k] - df.y[i]);
            int vi = f.c[i].v, vj = f.c[j].v, vk = f.c[k].v;
            double a = r[vi] + r[vj];
            double b = r[vi] + r[vk];
            double c = r[vj] + r[vk];
            double ideal = acos(clamp_unit((a * a + b * b - c * c) / (2.0 * a * b)));
            double err = fabs(actual - ideal);
            if (isfinite(err)) {
                dr.theta_vertex[vi] += err;
                dr.theta_vertex_count[vi]++;
                theta_sum += err;
                theta_sumsq += err * err;
                dr.theta_max_abs = isnan(dr.theta_max_abs) ? err : max(dr.theta_max_abs, err);
                theta_count++;
            }
        }
    }

    dr.global_vertex_max_spread = 0.0;
    for (const auto& kv : spread) dr.global_vertex_max_spread = max(dr.global_vertex_max_spread, kv.second.max_spread);
    dr.global_vertex_rms_spread = spread_count > 0 ? sqrt(spread_sumsq / spread_count) : 0.0;

    for (int i = 0; i < mesh.N; ++i) {
        if (dr.theta_vertex_count[i] > 0) dr.theta_vertex[i] /= static_cast<double>(dr.theta_vertex_count[i]);
    }
    if (theta_count > 0) {
        dr.theta_mean_abs = theta_sum / theta_count;
        dr.theta_rms = sqrt(theta_sumsq / theta_count);
        if (isnan(dr.theta_max_abs)) dr.theta_max_abs = 0.0;
    }

    double det = S00 * S11 - S01 * S01;
    if (constraints > 0 && fabs(det) > 1e-18) {
        dr.period_A_x = (rx0 * S11 - rx1 * S01) / det;
        dr.period_B_x = (S00 * rx1 - S01 * rx0) / det;
        dr.period_A_y = (ry0 * S11 - ry1 * S01) / det;
        dr.period_B_y = (S00 * ry1 - S01 * ry0) / det;

        double res_sumsq = 0.0;
        int res_count = 0;
        for (int fi = 0; fi < F; ++fi) {
            if (!dr.faces[fi].done) continue;
            const Face& f = mesh.faces[fi];
            const DevelopedFace& df = dr.faces[fi];
            for (int i = 0; i < 3; ++i) {
                int v = f.c[i].v;
                const RefOcc& ref = refs[v];
                double dsx = static_cast<double>(df.sx[i] - ref.sx);
                double dsy = static_cast<double>(df.sy[i] - ref.sy);
                if (dsx == 0.0 && dsy == 0.0) continue;
                double predx = dsx * dr.period_A_x + dsy * dr.period_B_x;
                double predy = dsx * dr.period_A_y + dsy * dr.period_B_y;
                double rx = (df.x[i] - ref.x) - predx;
                double ry = (df.y[i] - ref.y) - predy;
                res_sumsq += rx * rx + ry * ry;
                res_count++;
            }
        }
        dr.period_fit_rms = res_count > 0 ? sqrt(res_sumsq / res_count) : 0.0;
    }

    return dr;
}

static vector<int> parse_int_list(const string& s) {
    vector<int> out;
    string item;
    stringstream ss(s);
    while (getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(stoi(item));
    }
    return out;
}

static vector<unsigned long long> parse_seed_list(const string& s) {
    vector<unsigned long long> out;
    string item;
    stringstream ss(s);
    while (getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(stoull(item));
    }
    return out;
}

static bool parse_bool_int(const string& s) {
    if (s == "1" || s == "true" || s == "yes") return true;
    if (s == "0" || s == "false" || s == "no") return false;
    throw runtime_error("expected 0|1 boolean, got " + s);
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        auto val = [&](const string& name) -> string {
            if (i + 1 >= argc) throw runtime_error("missing value for " + name);
            return argv[++i];
        };
        if (arg == "--nx") cfg.nx = stoi(val(arg));
        else if (arg == "--ny") cfg.ny = stoi(val(arg));
        else if (arg == "--flips") cfg.flips = stoi(val(arg));
        else if (arg == "--seed") cfg.seed = stoull(val(arg));
        else if (arg == "--max_iter") cfg.max_iter = stoi(val(arg));
        else if (arg == "--tol") cfg.tol = stod(val(arg));
        else if (arg == "--method") cfg.method = val(arg);
        else if (arg == "--out") cfg.out = val(arg);
        else if (arg == "--bounded") cfg.bounded = parse_bool_int(val(arg));
        else if (arg == "--u_min") cfg.u_min = stod(val(arg));
        else if (arg == "--u_max") cfg.u_max = stod(val(arg));
        else if (arg == "--develop") cfg.develop = parse_bool_int(val(arg));
        else if (arg == "--embed") cfg.embed = parse_bool_int(val(arg));
        else if (arg == "--progress_every") cfg.progress_every = stoi(val(arg));
        else if (arg == "--validate_every") cfg.validate_every = stoi(val(arg));
        else if (arg == "--stagnation_window") cfg.stagnation_window = stoi(val(arg));
        else if (arg == "--stagnation_rel") cfg.stagnation_rel = stod(val(arg));
        else if (arg == "--scan_flips") cfg.scan_flips = parse_bool_int(val(arg));
        else if (arg == "--scan_list") cfg.scan_list = parse_int_list(val(arg));
        else if (arg == "--scan_seeds") cfg.scan_seeds = parse_seed_list(val(arg));
        else if (arg == "--scan_out") cfg.scan_out = val(arg);
        else if (arg == "--scan_write_details") cfg.scan_write_details = parse_bool_int(val(arg));
        else if (arg == "--threads") cfg.threads = stoi(val(arg));
        else if (arg == "--write_develop_faces") cfg.write_develop_faces = parse_bool_int(val(arg)) ? 1 : 0;
        else if (arg == "--write_theta_vertex") cfg.write_theta_vertex = parse_bool_int(val(arg)) ? 1 : 0;
        else if (arg == "--write_coords") cfg.write_coords = parse_bool_int(val(arg)) ? 1 : 0;
        else if (arg == "--test") cfg.test = true;
        else if (arg == "-h" || arg == "--help") {
            cout << "Usage: ./inc_ricci_v3 [--nx N] [--ny N] [--flips N] [--seed N] [--method ricci]\n";
            exit(0);
        } else {
            throw runtime_error("unknown argument: " + arg);
        }
    }
    if (cfg.nx <= 2 || cfg.ny <= 2) throw runtime_error("nx and ny must be > 2");
    if (cfg.flips < 0) throw runtime_error("flips must be >= 0");
    if (cfg.max_iter < 0) throw runtime_error("max_iter must be >= 0");
    if (cfg.tol <= 0.0) throw runtime_error("tol must be positive");
    if (cfg.method != "ricci") throw runtime_error("v3 production path supports --method ricci");
    if (cfg.u_min >= cfg.u_max) throw runtime_error("u_min must be < u_max");
    if (cfg.validate_every < 0) throw runtime_error("validate_every must be >= 0");
    if (cfg.progress_every < 0) throw runtime_error("progress_every must be >= 0");
    if (cfg.threads <= 0) cfg.threads = 1;
    return cfg;
}

static void apply_output_defaults(Config& cfg, bool scan_mode) {
    if (scan_mode && !cfg.scan_write_details) {
        if (cfg.write_develop_faces < 0) cfg.write_develop_faces = 0;
        if (cfg.write_theta_vertex < 0) cfg.write_theta_vertex = 0;
        if (cfg.write_coords < 0) cfg.write_coords = 0;
    } else {
        if (cfg.write_develop_faces < 0) cfg.write_develop_faces = 1;
        if (cfg.write_theta_vertex < 0) cfg.write_theta_vertex = 1;
        if (cfg.write_coords < 0) cfg.write_coords = 1;
    }
}

static map<int, int> degree_histogram(const Mesh& mesh) {
    vector<int> degree(mesh.N, 0);
    for (const EdgeAdj& e : mesh.edges) {
        degree[e.a]++;
        degree[e.b]++;
    }
    map<int, int> hist;
    for (int d : degree) hist[d]++;
    return hist;
}

static RunOutput execute_run(Config cfg, bool write_files, bool verbose) {
    apply_output_defaults(cfg, false);
    RunOutput out;
    out.cfg = cfg;
    double t_total0 = now_sec();
    try {
        double t0 = now_sec();
        mt19937_64 rng(cfg.seed);
        out.mesh = generate_initial_mesh(cfg.nx, cfg.ny);
        string error;
        if (!validate_mesh(out.mesh, &error)) throw runtime_error("initial mesh invalid: " + error);
        out.flips_accepted = randomize_topology(out.mesh, cfg.flips, cfg.validate_every, rng);
        if (!validate_mesh(out.mesh, &error)) throw runtime_error("final mesh invalid: " + error);
        out.time.topology = now_sec() - t0;

        vector<double> u = random_log_radii(out.mesh.N, rng, 0.1);
        t0 = now_sec();
        out.ricci = run_ricci(out.mesh, std::move(u), cfg, write_files, verbose);
        out.time.ricci = now_sec() - t0;

        if (cfg.develop) {
            t0 = now_sec();
            out.develop = develop_mesh(out.mesh, out.ricci.state.r);
            out.time.develop = now_sec() - t0;
        }

        out.time.total = now_sec() - t_total0;
    } catch (const exception& e) {
        out.valid = false;
        out.error = e.what();
        out.time.total = now_sec() - t_total0;
    }
    return out;
}

static void write_outputs(const RunOutput& out) {
    const Config& cfg = out.cfg;
    filesystem::create_directories(cfg.out);
    const Mesh& mesh = out.mesh;
    const RicciResult& rr = out.ricci;
    const DevelopResult& dr = out.develop;

    {
        ofstream f(filesystem::path(cfg.out) / "stats.csv");
        f << "iter,E_K,max_abs_K,step,r_min,r_max,accepted\n";
        f << scientific << setprecision(17);
        for (const StatRow& row : rr.stats) {
            f << row.iter << ',' << row.E_K << ',' << row.max_abs_K << ',' << row.step << ','
              << row.r_min << ',' << row.r_max << ',' << row.accepted << '\n';
        }
    }
    {
        ofstream f(filesystem::path(cfg.out) / "radii.csv");
        f << "vertex,u,r\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < mesh.N; ++i) f << i << ',' << rr.state.u[i] << ',' << rr.state.r[i] << '\n';
    }
    {
        ofstream f(filesystem::path(cfg.out) / "faces.dat");
        for (const Face& face : mesh.faces) f << face.c[0].v << ' ' << face.c[1].v << ' ' << face.c[2].v << '\n';
    }
    {
        ofstream f(filesystem::path(cfg.out) / "degree_hist.csv");
        f << "degree,count\n";
        for (const auto& kv : degree_histogram(mesh)) f << kv.first << ',' << kv.second << '\n';
    }
    {
        ofstream f(filesystem::path(cfg.out) / "edges.csv");
        f << "i,j,shift_x,shift_y,target_length\n";
        f << scientific << setprecision(17);
        for (const EdgeAdj& e : mesh.edges) {
            f << e.a << ',' << e.b << ',' << e.shift_x << ',' << e.shift_y << ','
              << rr.state.r[e.a] + rr.state.r[e.b] << '\n';
        }
    }
    if (cfg.write_develop_faces && dr.ran) {
        ofstream f(filesystem::path(cfg.out) / "develop_faces.csv");
        f << "face,corner,vertex,lift_sx,lift_sy,x,y\n";
        f << scientific << setprecision(17);
        for (int fi = 0; fi < static_cast<int>(mesh.faces.size()); ++fi) {
            if (!dr.faces[fi].done) continue;
            for (int q = 0; q < 3; ++q) {
                f << fi << ',' << q << ',' << mesh.faces[fi].c[q].v << ','
                  << dr.faces[fi].sx[q] << ',' << dr.faces[fi].sy[q] << ','
                  << dr.faces[fi].x[q] << ',' << dr.faces[fi].y[q] << '\n';
            }
        }
    }
    if (cfg.write_theta_vertex && dr.ran) {
        ofstream f(filesystem::path(cfg.out) / "theta_vertex.csv");
        f << "vertex,Theta_i,incident_face_corners\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < mesh.N; ++i) {
            f << i << ',' << dr.theta_vertex[i] << ',' << dr.theta_vertex_count[i] << '\n';
        }
    }
    if (cfg.write_coords && dr.ran) {
        vector<char> have(mesh.N, 0);
        vector<double> x(mesh.N, numeric_limits<double>::quiet_NaN());
        vector<double> y(mesh.N, numeric_limits<double>::quiet_NaN());
        for (int fi = 0; fi < static_cast<int>(mesh.faces.size()); ++fi) {
            if (!dr.faces[fi].done) continue;
            for (int q = 0; q < 3; ++q) {
                int v = mesh.faces[fi].c[q].v;
                if (!have[v]) {
                    have[v] = 1;
                    x[v] = dr.faces[fi].x[q];
                    y[v] = dr.faces[fi].y[q];
                }
            }
        }
        ofstream f(filesystem::path(cfg.out) / "coords.csv");
        f << "vertex,x_dev_ref,y_dev_ref,u,r\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < mesh.N; ++i) f << i << ',' << x[i] << ',' << y[i] << ',' << rr.state.u[i] << ',' << rr.state.r[i] << '\n';
    }
    if (dr.ran) {
        ofstream f(filesystem::path(cfg.out) / "develop_summary.txt");
        f << scientific << setprecision(17);
        f << "develop_status: " << dr.status << "\n";
        f << "develop_faces_count: " << dr.faces_count << "\n";
        f << "develop_global_vertex_max_spread: " << dr.global_vertex_max_spread << "\n";
        f << "develop_global_vertex_rms_spread: " << dr.global_vertex_rms_spread << "\n";
        f << "period_fit_rms: " << dr.period_fit_rms << "\n";
        f << "period_A: " << dr.period_A_x << " " << dr.period_A_y << "\n";
        f << "period_B: " << dr.period_B_x << " " << dr.period_B_y << "\n";
        f << "develop_theta_mean_abs: " << dr.theta_mean_abs << "\n";
        f << "develop_theta_rms: " << dr.theta_rms << "\n";
        f << "develop_theta_max_abs: " << dr.theta_max_abs << "\n";
    }
    {
        ofstream f(filesystem::path(cfg.out) / "summary.txt");
        f << scientific << setprecision(17);
        f << "INC-Ricci v3 performance-oriented combinatorial/development run\n";
        f << "N: " << mesh.N << "\n";
        f << "E: " << mesh.edges.size() << "\n";
        f << "F: " << mesh.faces.size() << "\n";
        f << "flips_requested: " << cfg.flips << "\n";
        f << "flips_accepted: " << out.flips_accepted << "\n";
        f << "seed: " << cfg.seed << "\n";
        f << "bounded: " << (cfg.bounded ? "yes" : "no") << "\n";
        f << "boundary_rejections: " << rr.boundary_rejections << "\n";
        f << "ricci_iters: " << rr.iterations << "\n";
        f << "final_E_K: " << rr.state.E_K << "\n";
        f << "final_max_abs_K: " << rr.state.max_abs_K << "\n";
        f << "radius_ratio: " << rr.radius_ratio << "\n";
        f << "stopped_reason: " << rr.stopped_reason << "\n";
        f << "develop_faces_count: " << dr.faces_count << "\n";
        f << "develop_global_vertex_max_spread: " << dr.global_vertex_max_spread << "\n";
        f << "period_fit_rms: " << dr.period_fit_rms << "\n";
        f << "develop_theta_mean_abs: " << dr.theta_mean_abs << "\n";
        f << "time_topology_sec: " << out.time.topology << "\n";
        f << "time_ricci_sec: " << out.time.ricci << "\n";
        f << "time_develop_sec: " << out.time.develop << "\n";
        f << "time_embed_sec: " << out.time.embed << "\n";
        f << "time_total_sec: " << out.time.total << "\n";
        f << "note: v3 development diagnostics are linear-time intrinsic unfolding diagnostics, not a full coordinate-level packing solver.\n";
    }
}

static ScanRow to_scan_row(const RunOutput& out) {
    ScanRow row;
    row.nx = out.cfg.nx;
    row.ny = out.cfg.ny;
    row.N = out.mesh.N;
    row.E = static_cast<int>(out.mesh.edges.size());
    row.F = static_cast<int>(out.mesh.faces.size());
    row.flips_requested = out.cfg.flips;
    row.flips_accepted = out.flips_accepted;
    row.seed = out.cfg.seed;
    row.ricci_iters = out.ricci.iterations;
    row.final_E_K = out.ricci.state.E_K;
    row.final_max_abs_K = out.ricci.state.max_abs_K;
    row.radius_ratio = out.ricci.radius_ratio;
    row.develop_global_vertex_max_spread = out.develop.global_vertex_max_spread;
    row.period_fit_rms = out.develop.period_fit_rms;
    row.develop_theta_mean_abs = out.develop.theta_mean_abs;
    row.time_total_sec = out.time.total;
    row.time_ricci_sec = out.time.ricci;
    row.time_develop_sec = out.time.develop;
    row.stopped_reason = out.ricci.stopped_reason;
    if (!out.valid) row.status = out.error;
    return row;
}

static void write_scan_csv(const string& path, const vector<ScanRow>& rows) {
    ofstream f(path);
    f << "nx,ny,N,E,F,flips_requested,flips_accepted,seed,ricci_iters,final_E_K,final_max_abs_K,"
         "radius_ratio,develop_global_vertex_max_spread,period_fit_rms,develop_theta_mean_abs,"
         "time_total_sec,time_ricci_sec,time_develop_sec,stopped_reason,status\n";
    f << scientific << setprecision(17);
    for (const ScanRow& r : rows) {
        f << r.nx << ',' << r.ny << ',' << r.N << ',' << r.E << ',' << r.F << ','
          << r.flips_requested << ',' << r.flips_accepted << ',' << r.seed << ','
          << r.ricci_iters << ',' << r.final_E_K << ',' << r.final_max_abs_K << ','
          << r.radius_ratio << ',' << r.develop_global_vertex_max_spread << ','
          << r.period_fit_rms << ',' << r.develop_theta_mean_abs << ','
          << r.time_total_sec << ',' << r.time_ricci_sec << ',' << r.time_develop_sec << ','
          << '"' << r.stopped_reason << '"' << ',' << '"' << r.status << '"' << '\n';
    }
}

static vector<ScanRow> run_scan(Config cfg) {
    apply_output_defaults(cfg, true);
    if (cfg.scan_list.empty()) cfg.scan_list.push_back(cfg.flips);
    if (cfg.scan_seeds.empty()) cfg.scan_seeds.push_back(cfg.seed);
    cfg.progress_every = 0;
    cfg.develop = true;

    struct Task { int flips; unsigned long long seed; };
    vector<Task> tasks;
    for (int flips : cfg.scan_list) {
        for (unsigned long long seed : cfg.scan_seeds) tasks.push_back({flips, seed});
    }

    vector<ScanRow> rows(tasks.size());
    atomic<size_t> next{0};
    int nthreads = max(1, min(cfg.threads, static_cast<int>(tasks.size())));
    vector<thread> workers;
    workers.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) {
        workers.emplace_back([&, t]() {
            while (true) {
                size_t idx = next.fetch_add(1);
                if (idx >= tasks.size()) break;
                Config local = cfg;
                local.flips = tasks[idx].flips;
                local.seed = tasks[idx].seed;
                RunOutput out = execute_run(local, cfg.scan_write_details, false);
                rows[idx] = to_scan_row(out);
            }
        });
    }
    for (thread& th : workers) th.join();

    for (const ScanRow& r : rows) {
        cout << scientific << setprecision(6)
             << "scan flips=" << r.flips_requested
             << " seed=" << r.seed
             << " maxK=" << r.final_max_abs_K
             << " spread=" << r.develop_global_vertex_max_spread
             << " time=" << r.time_total_sec
             << " status=" << r.status << '\n';
    }
    write_scan_csv(cfg.scan_out, rows);
    return rows;
}

static bool run_tests() {
    bool ok = true;
    cout << "Running v3 tests\n";
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 0; cfg.seed = 1; cfg.max_iter = 1;
        cfg.tol = 1e-12; cfg.develop = true; cfg.progress_every = 0;
        Mesh mesh = generate_initial_mesh(cfg.nx, cfg.ny);
        string error;
        bool valid = validate_mesh(mesh, &error);
        vector<double> u(mesh.N, 0.0);
        RicciResult rr = run_ricci(mesh, u, cfg, false, false);
        DevelopResult dr = develop_mesh(mesh, rr.state.r);
        bool pass = valid && rr.state.max_abs_K < 1e-12 && dr.faces_count == static_cast<int>(mesh.faces.size()) &&
                    dr.theta_mean_abs < 1e-12 && dr.period_fit_rms < 1e-10;
        ok = ok && pass;
        cout << "Test 1 regular flat develop: " << (pass ? "PASS" : "FAIL")
             << " maxK=" << scientific << setprecision(6) << rr.state.max_abs_K
             << " theta_mean=" << dr.theta_mean_abs
             << " period_fit=" << dr.period_fit_rms << '\n';
        if (!valid) cout << "  validation: " << error << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 200; cfg.seed = 1; cfg.max_iter = 10000;
        cfg.tol = 1e-8; cfg.develop = true; cfg.progress_every = 0;
        RunOutput out = execute_run(cfg, false, false);
        bool pass = out.valid && out.ricci.state.max_abs_K < 0.5 * out.ricci.initial_max_abs_K &&
                    out.ricci.state.E_K < out.ricci.initial_E_K && out.develop.faces_count == static_cast<int>(out.mesh.faces.size());
        ok = ok && pass;
        cout << "Test 2 random ricci/develop: " << (pass ? "PASS" : "FAIL")
             << " flips=" << out.flips_accepted
             << " initialK=" << out.ricci.initial_max_abs_K
             << " finalK=" << out.ricci.state.max_abs_K
             << " spread=" << out.develop.global_vertex_max_spread << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 200; cfg.seed = 1; cfg.max_iter = 10000;
        cfg.tol = 1e-8; cfg.bounded = true; cfg.u_min = -2.0; cfg.u_max = 2.0;
        cfg.develop = true; cfg.progress_every = 0;
        RunOutput out = execute_run(cfg, false, false);
        bool bounds_ok = out.ricci.u_min_actual >= cfg.u_min - 1e-12 && out.ricci.u_max_actual <= cfg.u_max + 1e-12;
        bool pass = out.valid && bounds_ok && out.ricci.state.E_K < out.ricci.initial_E_K;
        ok = ok && pass;
        cout << "Test 3 bounded ricci: " << (pass ? "PASS" : "FAIL")
             << " finalK=" << out.ricci.state.max_abs_K
             << " boundary_rejections=" << out.ricci.boundary_rejections
             << " u_range=[" << out.ricci.u_min_actual << "," << out.ricci.u_max_actual << "]\n";
    }
    cout << "Tests " << (ok ? "PASSED" : "FAILED") << '\n';
    return ok;
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);
        if (cfg.test) return run_tests() ? 0 : 1;

        if (cfg.scan_flips) {
            double t0 = now_sec();
            vector<ScanRow> rows = run_scan(cfg);
            double total = now_sec() - t0;
            cout << "scan_rows " << rows.size() << " scan_out " << cfg.scan_out
                 << " total_sec " << scientific << setprecision(6) << total
                 << " threads " << cfg.threads << '\n';
            return 0;
        }

        RunOutput out = execute_run(cfg, true, cfg.progress_every > 0);
        if (!out.valid) throw runtime_error(out.error);
        write_outputs(out);

        cout << scientific << setprecision(17);
        cout << "Done\n";
        cout << "out_dir " << cfg.out << '\n';
        cout << "N " << out.mesh.N << " E " << out.mesh.edges.size() << " F " << out.mesh.faces.size() << '\n';
        cout << "flips_requested " << cfg.flips << " flips_accepted " << out.flips_accepted << '\n';
        cout << "ricci_iters " << out.ricci.iterations << '\n';
        cout << "final_E_K " << out.ricci.state.E_K << '\n';
        cout << "final_max_abs_K " << out.ricci.state.max_abs_K << '\n';
        cout << "radius_ratio " << out.ricci.radius_ratio << '\n';
        cout << "develop_global_vertex_max_spread " << out.develop.global_vertex_max_spread << '\n';
        cout << "period_fit_rms " << out.develop.period_fit_rms << '\n';
        cout << "develop_theta_mean_abs " << out.develop.theta_mean_abs << '\n';
        cout << "time_topology_sec " << out.time.topology << '\n';
        cout << "time_ricci_sec " << out.time.ricci << '\n';
        cout << "time_develop_sec " << out.time.develop << '\n';
        cout << "time_total_sec " << out.time.total << '\n';
        cout << "stopped_reason " << out.ricci.stopped_reason << '\n';
        return 0;
    } catch (const exception& e) {
        cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
