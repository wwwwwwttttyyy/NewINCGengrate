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
#include <unordered_set>
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

struct RadiusStats {
    double mean_r = numeric_limits<double>::quiet_NaN();
    double std_r = numeric_limits<double>::quiet_NaN();
    double cv_r = numeric_limits<double>::quiet_NaN();
    double polydispersity_delta = numeric_limits<double>::quiet_NaN();
    double q01_r = numeric_limits<double>::quiet_NaN();
    double q05_r = numeric_limits<double>::quiet_NaN();
    double q50_r = numeric_limits<double>::quiet_NaN();
    double q95_r = numeric_limits<double>::quiet_NaN();
    double q99_r = numeric_limits<double>::quiet_NaN();
    double q01_u = numeric_limits<double>::quiet_NaN();
    double q05_u = numeric_limits<double>::quiet_NaN();
    double q50_u = numeric_limits<double>::quiet_NaN();
    double q95_u = numeric_limits<double>::quiet_NaN();
    double q99_u = numeric_limits<double>::quiet_NaN();
    vector<int> hist_counts;
    double hist_min = 0.0;
    double hist_max = 0.0;
};

struct OverlapSample {
    int i = 0;
    int j = 0;
    double dist = 0.0;
    double ri = 0.0;
    double rj = 0.0;
    double overlap = 0.0;
};

struct EdgeError {
    int i = 0;
    int j = 0;
    int shift_x = 0;
    int shift_y = 0;
    double dist = 0.0;
    double target = 0.0;
    double error = 0.0;
};

struct CrossingSample {
    int e1 = 0;
    int e2 = 0;
    int shift_x = 0;
    int shift_y = 0;
};

struct PackingDiagnostics {
    bool ran = false;
    bool cell_skew_warning = false;
    double cell_skew_sin = numeric_limits<double>::quiet_NaN();
    vector<double> base_x;
    vector<double> base_y;

    double max_edge_contact_error = numeric_limits<double>::quiet_NaN();
    double rms_edge_contact_error = numeric_limits<double>::quiet_NaN();
    int nonedge_overlap_count = 0;
    double nonedge_max_overlap = numeric_limits<double>::quiet_NaN();
    double nonedge_rms_overlap = numeric_limits<double>::quiet_NaN();
    double nonedge_overlap_fraction = numeric_limits<double>::quiet_NaN();
    vector<OverlapSample> overlap_samples;
    vector<EdgeError> edge_errors;

    double face_min_area = numeric_limits<double>::quiet_NaN();
    double face_max_area = numeric_limits<double>::quiet_NaN();
    double face_mean_area = numeric_limits<double>::quiet_NaN();
    int face_bad_area_count = 0;

    long long edge_crossing_count = 0;
    long long edge_crossing_checked_pairs = 0;
    double edge_crossing_fraction = numeric_limits<double>::quiet_NaN();
    vector<CrossingSample> crossing_samples;

    RadiusStats radius_stats;

    bool intrinsic_endpoint = false;
    bool developed_consistent = false;
    bool local_theta_valid = false;
    bool contact_edges_valid = false;
    bool nonedge_overlap_free = false;
    bool topology_geometry_valid = false;
    bool physically_reasonable_radius = false;
    bool valid_disk_packing_candidate = false;
    bool strong_INC_candidate = false;
    bool physical_INC_like_candidate = false;
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
    bool scan_bounds = false;
    vector<int> scan_list;
    vector<unsigned long long> scan_seeds;
    string scan_out = "scan.csv";
    string bounds_out = "bounds_scan.csv";
    vector<double> bounds_list;
    bool include_unbounded = false;
    bool scan_write_details = false;
    int threads = 4;

    int write_develop_faces = -1;
    int write_theta_vertex = -1;
    int write_coords = -1;

    double overlap_tol = 1e-8;
    string overlap_mode = "full";
    double area_tol = 1e-12;
    bool check_crossings = true;
    bool check_crossings_explicit = false;
    double intersect_tol = 1e-10;
    int hist_bins = 50;
    double target_max_ratio = 10.0;
    double target_delta = 0.25;
    double develop_spread_tol = 1e-6;
    double period_fit_tol = 1e-6;
    double theta_tol = 1e-8;
    double contact_tol = 1e-6;

    bool surgery = false;
    int surgery_steps = 1000;
    int surgery_trials_per_step = 20;
    string surgery_metric = "overlap";
    double surgery_temperature = 0.0;
    int surgery_ricci_iter = 20000;
    int surgery_report_every = 10;
};

struct RunOutput {
    Config cfg;
    Mesh mesh;
    RicciResult ricci;
    DevelopResult develop;
    PackingDiagnostics packing;
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
    int nonedge_overlap_count = 0;
    double nonedge_max_overlap = numeric_limits<double>::quiet_NaN();
    long long edge_crossing_count = 0;
    int face_bad_area_count = 0;
    double polydispersity_delta = numeric_limits<double>::quiet_NaN();
    bool physically_reasonable_radius = false;
    bool valid_disk_packing_candidate = false;
    bool strong_INC_candidate = false;
    double time_total_sec = 0.0;
    double time_ricci_sec = 0.0;
    double time_develop_sec = 0.0;
    string stopped_reason;
    string status = "ok";
};

struct BoundsRow {
    int nx = 0, ny = 0;
    unsigned long long seed = 0;
    int flips_requested = 0;
    int flips_accepted = 0;
    double B = numeric_limits<double>::quiet_NaN();
    int unbounded = 0;
    double final_E_K = numeric_limits<double>::quiet_NaN();
    double final_max_abs_K = numeric_limits<double>::quiet_NaN();
    double radius_ratio = numeric_limits<double>::quiet_NaN();
    double polydispersity_delta = numeric_limits<double>::quiet_NaN();
    int boundary_rejections = 0;
    int ricci_iters = 0;
    double develop_global_vertex_max_spread = numeric_limits<double>::quiet_NaN();
    double period_fit_rms = numeric_limits<double>::quiet_NaN();
    int nonedge_overlap_count = 0;
    double nonedge_max_overlap = numeric_limits<double>::quiet_NaN();
    long long edge_crossing_count = 0;
    int face_bad_area_count = 0;
    bool physically_reasonable_radius = false;
    bool valid_disk_packing_candidate = false;
    bool strong_INC_candidate = false;
    double time_total_sec = 0.0;
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

static vector<double> parse_double_list(const string& s) {
    vector<double> out;
    string item;
    stringstream ss(s);
    while (getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(stod(item));
    }
    return out;
}

static double quantile_sorted(const vector<double>& x, double q) {
    if (x.empty()) return numeric_limits<double>::quiet_NaN();
    double pos = q * static_cast<double>(x.size() - 1);
    size_t lo = static_cast<size_t>(floor(pos));
    size_t hi = min(lo + 1, x.size() - 1);
    double t = pos - static_cast<double>(lo);
    return (1.0 - t) * x[lo] + t * x[hi];
}

static RadiusStats compute_radius_stats(const vector<double>& r, const vector<double>& u, int bins) {
    RadiusStats rs;
    if (r.empty()) return rs;
    double sum = accumulate(r.begin(), r.end(), 0.0);
    double sum2 = 0.0;
    for (double x : r) sum2 += x * x;
    rs.mean_r = sum / static_cast<double>(r.size());
    double var = max(0.0, sum2 / static_cast<double>(r.size()) - rs.mean_r * rs.mean_r);
    rs.std_r = sqrt(var);
    rs.cv_r = rs.mean_r > 0.0 ? rs.std_r / rs.mean_r : numeric_limits<double>::infinity();
    rs.polydispersity_delta = rs.cv_r;

    vector<double> sr = r;
    vector<double> su = u;
    sort(sr.begin(), sr.end());
    sort(su.begin(), su.end());
    rs.q01_r = quantile_sorted(sr, 0.01);
    rs.q05_r = quantile_sorted(sr, 0.05);
    rs.q50_r = quantile_sorted(sr, 0.50);
    rs.q95_r = quantile_sorted(sr, 0.95);
    rs.q99_r = quantile_sorted(sr, 0.99);
    rs.q01_u = quantile_sorted(su, 0.01);
    rs.q05_u = quantile_sorted(su, 0.05);
    rs.q50_u = quantile_sorted(su, 0.50);
    rs.q95_u = quantile_sorted(su, 0.95);
    rs.q99_u = quantile_sorted(su, 0.99);

    bins = max(1, bins);
    rs.hist_counts.assign(bins, 0);
    rs.hist_min = sr.front();
    rs.hist_max = sr.back();
    double width = rs.hist_max - rs.hist_min;
    for (double x : r) {
        int b = width > 0.0 ? static_cast<int>((x - rs.hist_min) / width * bins) : 0;
        if (b < 0) b = 0;
        if (b >= bins) b = bins - 1;
        rs.hist_counts[b]++;
    }
    return rs;
}

static bool reconstruct_base_positions(const Mesh& mesh, const DevelopResult& dr,
                                       vector<double>& px, vector<double>& py) {
    px.assign(mesh.N, 0.0);
    py.assign(mesh.N, 0.0);
    vector<int> count(mesh.N, 0);
    if (!dr.ran || !isfinite(dr.period_A_x) || !isfinite(dr.period_B_x)) return false;
    for (int fi = 0; fi < static_cast<int>(mesh.faces.size()); ++fi) {
        if (!dr.faces[fi].done) continue;
        const Face& f = mesh.faces[fi];
        const DevelopedFace& df = dr.faces[fi];
        for (int q = 0; q < 3; ++q) {
            int v = f.c[q].v;
            double bx = df.x[q] - static_cast<double>(df.sx[q]) * dr.period_A_x
                       - static_cast<double>(df.sy[q]) * dr.period_B_x;
            double by = df.y[q] - static_cast<double>(df.sx[q]) * dr.period_A_y
                       - static_cast<double>(df.sy[q]) * dr.period_B_y;
            px[v] += bx;
            py[v] += by;
            count[v]++;
        }
    }
    bool ok = true;
    for (int i = 0; i < mesh.N; ++i) {
        if (count[i] == 0) {
            ok = false;
        } else {
            px[i] /= static_cast<double>(count[i]);
            py[i] /= static_cast<double>(count[i]);
        }
    }
    return ok;
}

static double periodic_min_dist(double xi, double yi, double xj, double yj,
                                double Ax, double Ay, double Bx, double By) {
    double best2 = numeric_limits<double>::infinity();
    for (int m = -1; m <= 1; ++m) {
        for (int n = -1; n <= 1; ++n) {
            double dx = xj + m * Ax + n * Bx - xi;
            double dy = yj + m * Ay + n * By - yi;
            best2 = min(best2, dx * dx + dy * dy);
        }
    }
    return sqrt(best2);
}

static double orient2d(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

static bool on_segment_tol(double ax, double ay, double bx, double by, double px, double py, double tol) {
    return fabs(orient2d(ax, ay, bx, by, px, py)) <= tol &&
           px >= min(ax, bx) - tol && px <= max(ax, bx) + tol &&
           py >= min(ay, by) - tol && py <= max(ay, by) + tol;
}

static bool segments_intersect_tol(double ax, double ay, double bx, double by,
                                   double cx, double cy, double dx, double dy, double tol) {
    double o1 = orient2d(ax, ay, bx, by, cx, cy);
    double o2 = orient2d(ax, ay, bx, by, dx, dy);
    double o3 = orient2d(cx, cy, dx, dy, ax, ay);
    double o4 = orient2d(cx, cy, dx, dy, bx, by);
    if (((o1 > tol && o2 < -tol) || (o1 < -tol && o2 > tol)) &&
        ((o3 > tol && o4 < -tol) || (o3 < -tol && o4 > tol))) {
        return true;
    }
    if (on_segment_tol(ax, ay, bx, by, cx, cy, tol)) return true;
    if (on_segment_tol(ax, ay, bx, by, dx, dy, tol)) return true;
    if (on_segment_tol(cx, cy, dx, dy, ax, ay, tol)) return true;
    if (on_segment_tol(cx, cy, dx, dy, bx, by, tol)) return true;
    return false;
}

static PackingDiagnostics compute_packing_diagnostics(const Mesh& mesh, const RicciResult& rr,
                                                      const DevelopResult& dr, const Config& cfg) {
    PackingDiagnostics pd;
    pd.ran = true;
    pd.radius_stats = compute_radius_stats(rr.state.r, rr.state.u, cfg.hist_bins);
    pd.intrinsic_endpoint = rr.state.max_abs_K < cfg.tol;
    pd.developed_consistent = dr.ran &&
                              dr.global_vertex_max_spread < cfg.develop_spread_tol &&
                              dr.period_fit_rms < cfg.period_fit_tol;
    pd.local_theta_valid = dr.ran && dr.theta_mean_abs < cfg.theta_tol;
    pd.physically_reasonable_radius = rr.radius_ratio <= cfg.target_max_ratio &&
                                      pd.radius_stats.polydispersity_delta <= cfg.target_delta;

    if (!dr.ran || !reconstruct_base_positions(mesh, dr, pd.base_x, pd.base_y)) {
        pd.contact_edges_valid = false;
        pd.nonedge_overlap_free = false;
        pd.topology_geometry_valid = false;
        return pd;
    }

    double Ax = dr.period_A_x, Ay = dr.period_A_y, Bx = dr.period_B_x, By = dr.period_B_y;
    double area_cell = fabs(Ax * By - Ay * Bx);
    double lenA = sqrt(Ax * Ax + Ay * Ay);
    double lenB = sqrt(Bx * Bx + By * By);
    pd.cell_skew_sin = (lenA > 0.0 && lenB > 0.0) ? area_cell / (lenA * lenB) : 0.0;
    pd.cell_skew_warning = pd.cell_skew_sin < 0.1;

    unordered_set<uint64_t> edge_set;
    edge_set.reserve(mesh.edges.size() * 2 + 1);
    double edge_sumsq = 0.0;
    pd.max_edge_contact_error = 0.0;
    pd.edge_errors.reserve(mesh.edges.size());
    for (const EdgeAdj& e : mesh.edges) {
        edge_set.insert(pack_edge_key(e.a, e.b));
        double xj = pd.base_x[e.b] + e.shift_x * Ax + e.shift_y * Bx;
        double yj = pd.base_y[e.b] + e.shift_x * Ay + e.shift_y * By;
        double dx = xj - pd.base_x[e.a];
        double dy = yj - pd.base_y[e.a];
        double dist = sqrt(dx * dx + dy * dy);
        double target = rr.state.r[e.a] + rr.state.r[e.b];
        double err = dist - target;
        pd.max_edge_contact_error = max(pd.max_edge_contact_error, fabs(err));
        edge_sumsq += err * err;
        pd.edge_errors.push_back({e.a, e.b, e.shift_x, e.shift_y, dist, target, err});
    }
    pd.rms_edge_contact_error = mesh.edges.empty() ? 0.0 : sqrt(edge_sumsq / mesh.edges.size());
    pd.contact_edges_valid = pd.max_edge_contact_error < cfg.contact_tol;

    long long nonedge_pairs = 0;
    double overlap_sumsq = 0.0;
    pd.nonedge_max_overlap = 0.0;
    for (int i = 0; i < mesh.N; ++i) {
        for (int j = i + 1; j < mesh.N; ++j) {
            if (edge_set.find(pack_edge_key(i, j)) != edge_set.end()) continue;
            nonedge_pairs++;
            double dist = periodic_min_dist(pd.base_x[i], pd.base_y[i], pd.base_x[j], pd.base_y[j], Ax, Ay, Bx, By);
            double overlap = max(0.0, rr.state.r[i] + rr.state.r[j] - dist);
            if (overlap > cfg.overlap_tol) pd.nonedge_overlap_count++;
            if (overlap > 0.0) overlap_sumsq += overlap * overlap;
            pd.nonedge_max_overlap = max(pd.nonedge_max_overlap, overlap);
            if (overlap > 0.0) pd.overlap_samples.push_back({i, j, dist, rr.state.r[i], rr.state.r[j], overlap});
        }
    }
    pd.nonedge_overlap_fraction = nonedge_pairs > 0 ? static_cast<double>(pd.nonedge_overlap_count) / nonedge_pairs : 0.0;
    pd.nonedge_rms_overlap = nonedge_pairs > 0 ? sqrt(overlap_sumsq / nonedge_pairs) : 0.0;
    sort(pd.overlap_samples.begin(), pd.overlap_samples.end(),
         [](const OverlapSample& a, const OverlapSample& b) { return a.overlap > b.overlap; });
    if (pd.overlap_samples.size() > 100) pd.overlap_samples.resize(100);
    pd.nonedge_overlap_free = pd.nonedge_max_overlap < cfg.overlap_tol;

    pd.face_min_area = numeric_limits<double>::infinity();
    pd.face_max_area = -numeric_limits<double>::infinity();
    double area_sum = 0.0;
    int area_count = 0;
    for (int fi = 0; fi < static_cast<int>(mesh.faces.size()); ++fi) {
        if (!dr.faces[fi].done) continue;
        const DevelopedFace& df = dr.faces[fi];
        double area = 0.5 * orient2d(df.x[0], df.y[0], df.x[1], df.y[1], df.x[2], df.y[2]);
        if (area <= cfg.area_tol) pd.face_bad_area_count++;
        pd.face_min_area = min(pd.face_min_area, area);
        pd.face_max_area = max(pd.face_max_area, area);
        area_sum += area;
        area_count++;
    }
    if (area_count > 0) {
        pd.face_mean_area = area_sum / area_count;
    } else {
        pd.face_min_area = pd.face_max_area = pd.face_mean_area = numeric_limits<double>::quiet_NaN();
    }

    if (cfg.check_crossings) {
        for (int e1 = 0; e1 < static_cast<int>(mesh.edges.size()); ++e1) {
            const EdgeAdj& a = mesh.edges[e1];
            double a0x = pd.base_x[a.a], a0y = pd.base_y[a.a];
            double a1x = pd.base_x[a.b] + a.shift_x * Ax + a.shift_y * Bx;
            double a1y = pd.base_y[a.b] + a.shift_x * Ay + a.shift_y * By;
            for (int e2 = e1 + 1; e2 < static_cast<int>(mesh.edges.size()); ++e2) {
                const EdgeAdj& b = mesh.edges[e2];
                if (a.a == b.a || a.a == b.b || a.b == b.a || a.b == b.b) continue;
                for (int m = -1; m <= 1; ++m) {
                    for (int n = -1; n <= 1; ++n) {
                        double tx = m * Ax + n * Bx;
                        double ty = m * Ay + n * By;
                        double b0x = pd.base_x[b.a] + tx;
                        double b0y = pd.base_y[b.a] + ty;
                        double b1x = pd.base_x[b.b] + b.shift_x * Ax + b.shift_y * Bx + tx;
                        double b1y = pd.base_y[b.b] + b.shift_x * Ay + b.shift_y * By + ty;
                        pd.edge_crossing_checked_pairs++;
                        if (segments_intersect_tol(a0x, a0y, a1x, a1y, b0x, b0y, b1x, b1y, cfg.intersect_tol)) {
                            pd.edge_crossing_count++;
                            if (pd.crossing_samples.size() < 50) pd.crossing_samples.push_back({e1, e2, m, n});
                        }
                    }
                }
            }
        }
        pd.edge_crossing_fraction = pd.edge_crossing_checked_pairs > 0
                                      ? static_cast<double>(pd.edge_crossing_count) / pd.edge_crossing_checked_pairs
                                      : 0.0;
    } else {
        pd.edge_crossing_count = 0;
        pd.edge_crossing_checked_pairs = 0;
        pd.edge_crossing_fraction = numeric_limits<double>::quiet_NaN();
    }

    pd.topology_geometry_valid = pd.face_bad_area_count == 0 &&
                                 (!cfg.check_crossings || pd.edge_crossing_count == 0);
    pd.valid_disk_packing_candidate = pd.intrinsic_endpoint &&
                                      pd.developed_consistent &&
                                      pd.contact_edges_valid &&
                                      pd.nonedge_overlap_free &&
                                      pd.topology_geometry_valid;
    pd.strong_INC_candidate = pd.valid_disk_packing_candidate && pd.local_theta_valid;
    pd.physical_INC_like_candidate = pd.strong_INC_candidate && pd.physically_reasonable_radius;
    return pd;
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
        else if (arg == "--scan_bounds") cfg.scan_bounds = parse_bool_int(val(arg));
        else if (arg == "--scan_list") cfg.scan_list = parse_int_list(val(arg));
        else if (arg == "--scan_seeds") cfg.scan_seeds = parse_seed_list(val(arg));
        else if (arg == "--scan_out") cfg.scan_out = val(arg);
        else if (arg == "--bounds_out") cfg.bounds_out = val(arg);
        else if (arg == "--bounds_list") cfg.bounds_list = parse_double_list(val(arg));
        else if (arg == "--include_unbounded") cfg.include_unbounded = parse_bool_int(val(arg));
        else if (arg == "--scan_write_details") cfg.scan_write_details = parse_bool_int(val(arg));
        else if (arg == "--threads") cfg.threads = stoi(val(arg));
        else if (arg == "--write_develop_faces") cfg.write_develop_faces = parse_bool_int(val(arg)) ? 1 : 0;
        else if (arg == "--write_theta_vertex") cfg.write_theta_vertex = parse_bool_int(val(arg)) ? 1 : 0;
        else if (arg == "--write_coords") cfg.write_coords = parse_bool_int(val(arg)) ? 1 : 0;
        else if (arg == "--overlap_tol") cfg.overlap_tol = stod(val(arg));
        else if (arg == "--overlap_mode") cfg.overlap_mode = val(arg);
        else if (arg == "--area_tol") cfg.area_tol = stod(val(arg));
        else if (arg == "--check_crossings") {
            cfg.check_crossings = parse_bool_int(val(arg));
            cfg.check_crossings_explicit = true;
        }
        else if (arg == "--intersect_tol") cfg.intersect_tol = stod(val(arg));
        else if (arg == "--hist_bins") cfg.hist_bins = stoi(val(arg));
        else if (arg == "--target_max_ratio") cfg.target_max_ratio = stod(val(arg));
        else if (arg == "--target_delta") cfg.target_delta = stod(val(arg));
        else if (arg == "--develop_spread_tol") cfg.develop_spread_tol = stod(val(arg));
        else if (arg == "--period_fit_tol") cfg.period_fit_tol = stod(val(arg));
        else if (arg == "--theta_tol") cfg.theta_tol = stod(val(arg));
        else if (arg == "--contact_tol") cfg.contact_tol = stod(val(arg));
        else if (arg == "--surgery") cfg.surgery = parse_bool_int(val(arg));
        else if (arg == "--surgery_steps") cfg.surgery_steps = stoi(val(arg));
        else if (arg == "--surgery_trials_per_step") cfg.surgery_trials_per_step = stoi(val(arg));
        else if (arg == "--surgery_metric") cfg.surgery_metric = val(arg);
        else if (arg == "--surgery_temperature") cfg.surgery_temperature = stod(val(arg));
        else if (arg == "--surgery_ricci_iter") cfg.surgery_ricci_iter = stoi(val(arg));
        else if (arg == "--surgery_report_every") cfg.surgery_report_every = stoi(val(arg));
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
    if (cfg.hist_bins <= 0) throw runtime_error("hist_bins must be > 0");
    if (cfg.overlap_mode != "full" && cfg.overlap_mode != "cell") throw runtime_error("overlap_mode must be full or cell");
    if (cfg.surgery_metric != "overlap" && cfg.surgery_metric != "radius" && cfg.surgery_metric != "mixed") {
        throw runtime_error("surgery_metric must be overlap, radius, or mixed");
    }
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
            out.packing = compute_packing_diagnostics(out.mesh, out.ricci, out.develop, cfg);
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
    const PackingDiagnostics& pd = out.packing;

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
    if (pd.ran) {
        {
            ofstream f(filesystem::path(cfg.out) / "radius_stats.txt");
            f << scientific << setprecision(17);
            f << "r_min: " << rr.state.r_min << "\n";
            f << "r_max: " << rr.state.r_max << "\n";
            f << "radius_ratio: " << rr.radius_ratio << "\n";
            f << "mean_r: " << pd.radius_stats.mean_r << "\n";
            f << "std_r: " << pd.radius_stats.std_r << "\n";
            f << "cv_r: " << pd.radius_stats.cv_r << "\n";
            f << "polydispersity_delta: " << pd.radius_stats.polydispersity_delta << "\n";
            f << "q01_r: " << pd.radius_stats.q01_r << "\n";
            f << "q05_r: " << pd.radius_stats.q05_r << "\n";
            f << "q50_r: " << pd.radius_stats.q50_r << "\n";
            f << "q95_r: " << pd.radius_stats.q95_r << "\n";
            f << "q99_r: " << pd.radius_stats.q99_r << "\n";
            f << "q01_u: " << pd.radius_stats.q01_u << "\n";
            f << "q05_u: " << pd.radius_stats.q05_u << "\n";
            f << "q50_u: " << pd.radius_stats.q50_u << "\n";
            f << "q95_u: " << pd.radius_stats.q95_u << "\n";
            f << "q99_u: " << pd.radius_stats.q99_u << "\n";
        }
        {
            ofstream f(filesystem::path(cfg.out) / "radius_hist.csv");
            f << "bin,r_low,r_high,count\n";
            f << scientific << setprecision(17);
            int bins = static_cast<int>(pd.radius_stats.hist_counts.size());
            double width = bins > 0 ? (pd.radius_stats.hist_max - pd.radius_stats.hist_min) / bins : 0.0;
            for (int b = 0; b < bins; ++b) {
                double lo = pd.radius_stats.hist_min + b * width;
                double hi = (b == bins - 1) ? pd.radius_stats.hist_max : lo + width;
                f << b << ',' << lo << ',' << hi << ',' << pd.radius_stats.hist_counts[b] << '\n';
            }
        }
        {
            ofstream f(filesystem::path(cfg.out) / "packing_diagnostics.txt");
            f << scientific << setprecision(17);
            f << "intrinsic_endpoint: " << (pd.intrinsic_endpoint ? "yes" : "no") << "\n";
            f << "developed_consistent: " << (pd.developed_consistent ? "yes" : "no") << "\n";
            f << "local_theta_valid: " << (pd.local_theta_valid ? "yes" : "no") << "\n";
            f << "contact_edges_valid: " << (pd.contact_edges_valid ? "yes" : "no") << "\n";
            f << "nonedge_overlap_free: " << (pd.nonedge_overlap_free ? "yes" : "no") << "\n";
            f << "topology_geometry_valid: " << (pd.topology_geometry_valid ? "yes" : "no") << "\n";
            f << "physically_reasonable_radius: " << (pd.physically_reasonable_radius ? "yes" : "no") << "\n";
            f << "valid_disk_packing_candidate: " << (pd.valid_disk_packing_candidate ? "yes" : "no") << "\n";
            f << "strong_INC_candidate: " << (pd.strong_INC_candidate ? "yes" : "no") << "\n";
            f << "physical_INC_like_candidate: " << (pd.physical_INC_like_candidate ? "yes" : "no") << "\n";
            f << "max_edge_contact_error: " << pd.max_edge_contact_error << "\n";
            f << "rms_edge_contact_error: " << pd.rms_edge_contact_error << "\n";
            f << "nonedge_overlap_count: " << pd.nonedge_overlap_count << "\n";
            f << "nonedge_max_overlap: " << pd.nonedge_max_overlap << "\n";
            f << "nonedge_rms_overlap: " << pd.nonedge_rms_overlap << "\n";
            f << "nonedge_overlap_fraction: " << pd.nonedge_overlap_fraction << "\n";
            f << "face_min_area: " << pd.face_min_area << "\n";
            f << "face_max_area: " << pd.face_max_area << "\n";
            f << "face_mean_area: " << pd.face_mean_area << "\n";
            f << "face_bad_area_count: " << pd.face_bad_area_count << "\n";
            f << "edge_crossing_count: " << pd.edge_crossing_count << "\n";
            f << "edge_crossing_fraction: " << pd.edge_crossing_fraction << "\n";
            f << "cell_skew_sin: " << pd.cell_skew_sin << "\n";
            f << "cell_skew_warning: " << (pd.cell_skew_warning ? "yes" : "no") << "\n";
        }
        {
            ofstream f(filesystem::path(cfg.out) / "nonedge_overlaps_sample.csv");
            f << "i,j,dist,ri,rj,overlap_depth\n";
            f << scientific << setprecision(17);
            for (const OverlapSample& s : pd.overlap_samples) {
                f << s.i << ',' << s.j << ',' << s.dist << ',' << s.ri << ',' << s.rj << ',' << s.overlap << '\n';
            }
        }
        {
            ofstream f(filesystem::path(cfg.out) / "edge_contact_errors.csv");
            f << "i,j,shift_x,shift_y,dist,target,error\n";
            f << scientific << setprecision(17);
            for (const EdgeError& e : pd.edge_errors) {
                f << e.i << ',' << e.j << ',' << e.shift_x << ',' << e.shift_y << ','
                  << e.dist << ',' << e.target << ',' << e.error << '\n';
            }
        }
        {
            ofstream f(filesystem::path(cfg.out) / "crossings_sample.txt");
            for (const CrossingSample& c : pd.crossing_samples) {
                f << "edge1=" << c.e1 << " edge2=" << c.e2
                  << " shift=(" << c.shift_x << "," << c.shift_y << ")\n";
            }
        }
    }
    {
        ofstream f(filesystem::path(cfg.out) / "summary.txt");
        f << scientific << setprecision(17);
        f << "INC-Ricci v4 packing-diagnostic combinatorial/development run\n";
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
        f << "polydispersity_delta: " << pd.radius_stats.polydispersity_delta << "\n";
        f << "max_edge_contact_error: " << pd.max_edge_contact_error << "\n";
        f << "nonedge_overlap_count: " << pd.nonedge_overlap_count << "\n";
        f << "nonedge_max_overlap: " << pd.nonedge_max_overlap << "\n";
        f << "edge_crossing_count: " << pd.edge_crossing_count << "\n";
        f << "face_bad_area_count: " << pd.face_bad_area_count << "\n";
        f << "valid_disk_packing_candidate: " << (pd.valid_disk_packing_candidate ? "yes" : "no") << "\n";
        f << "strong_INC_candidate: " << (pd.strong_INC_candidate ? "yes" : "no") << "\n";
        f << "physical_INC_like_candidate: " << (pd.physical_INC_like_candidate ? "yes" : "no") << "\n";
        f << "time_topology_sec: " << out.time.topology << "\n";
        f << "time_ricci_sec: " << out.time.ricci << "\n";
        f << "time_develop_sec: " << out.time.develop << "\n";
        f << "time_embed_sec: " << out.time.embed << "\n";
        f << "time_total_sec: " << out.time.total << "\n";
        f << "note: v4 diagnostics test developed disk-packing geometry, but they are not a full coordinate-level packing solver.\n";
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
    row.nonedge_overlap_count = out.packing.nonedge_overlap_count;
    row.nonedge_max_overlap = out.packing.nonedge_max_overlap;
    row.edge_crossing_count = out.packing.edge_crossing_count;
    row.face_bad_area_count = out.packing.face_bad_area_count;
    row.polydispersity_delta = out.packing.radius_stats.polydispersity_delta;
    row.physically_reasonable_radius = out.packing.physically_reasonable_radius;
    row.valid_disk_packing_candidate = out.packing.valid_disk_packing_candidate;
    row.strong_INC_candidate = out.packing.strong_INC_candidate;
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
         "radius_ratio,polydispersity_delta,develop_global_vertex_max_spread,period_fit_rms,develop_theta_mean_abs,"
         "nonedge_overlap_count,nonedge_max_overlap,edge_crossing_count,face_bad_area_count,"
         "physically_reasonable_radius,valid_disk_packing_candidate,strong_INC_candidate,"
         "time_total_sec,time_ricci_sec,time_develop_sec,stopped_reason,status\n";
    f << scientific << setprecision(17);
    for (const ScanRow& r : rows) {
        f << r.nx << ',' << r.ny << ',' << r.N << ',' << r.E << ',' << r.F << ','
          << r.flips_requested << ',' << r.flips_accepted << ',' << r.seed << ','
          << r.ricci_iters << ',' << r.final_E_K << ',' << r.final_max_abs_K << ','
          << r.radius_ratio << ',' << r.polydispersity_delta << ',' << r.develop_global_vertex_max_spread << ','
          << r.period_fit_rms << ',' << r.develop_theta_mean_abs << ','
          << r.nonedge_overlap_count << ',' << r.nonedge_max_overlap << ',' << r.edge_crossing_count << ','
          << r.face_bad_area_count << ','
          << (r.physically_reasonable_radius ? 1 : 0) << ','
          << (r.valid_disk_packing_candidate ? 1 : 0) << ','
          << (r.strong_INC_candidate ? 1 : 0) << ','
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
    if (!cfg.check_crossings_explicit) cfg.check_crossings = false;

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

static BoundsRow to_bounds_row(const RunOutput& out, double B, int unbounded) {
    BoundsRow r;
    r.nx = out.cfg.nx;
    r.ny = out.cfg.ny;
    r.seed = out.cfg.seed;
    r.flips_requested = out.cfg.flips;
    r.flips_accepted = out.flips_accepted;
    r.B = B;
    r.unbounded = unbounded;
    r.final_E_K = out.ricci.state.E_K;
    r.final_max_abs_K = out.ricci.state.max_abs_K;
    r.radius_ratio = out.ricci.radius_ratio;
    r.polydispersity_delta = out.packing.radius_stats.polydispersity_delta;
    r.boundary_rejections = out.ricci.boundary_rejections;
    r.ricci_iters = out.ricci.iterations;
    r.develop_global_vertex_max_spread = out.develop.global_vertex_max_spread;
    r.period_fit_rms = out.develop.period_fit_rms;
    r.nonedge_overlap_count = out.packing.nonedge_overlap_count;
    r.nonedge_max_overlap = out.packing.nonedge_max_overlap;
    r.edge_crossing_count = out.packing.edge_crossing_count;
    r.face_bad_area_count = out.packing.face_bad_area_count;
    r.physically_reasonable_radius = out.packing.physically_reasonable_radius;
    r.valid_disk_packing_candidate = out.packing.valid_disk_packing_candidate;
    r.strong_INC_candidate = out.packing.strong_INC_candidate;
    r.time_total_sec = out.time.total;
    if (!out.valid) r.status = out.error;
    return r;
}

static void write_bounds_csv(const string& path, const vector<BoundsRow>& rows) {
    ofstream f(path);
    f << "nx,ny,seed,flips_requested,flips_accepted,B,unbounded,final_E_K,final_max_abs_K,"
         "radius_ratio,polydispersity_delta,boundary_rejections,ricci_iters,"
         "develop_global_vertex_max_spread,period_fit_rms,nonedge_overlap_count,nonedge_max_overlap,"
         "edge_crossing_count,face_bad_area_count,physically_reasonable_radius,valid_disk_packing_candidate,"
         "strong_INC_candidate,time_total_sec,status\n";
    f << scientific << setprecision(17);
    for (const BoundsRow& r : rows) {
        f << r.nx << ',' << r.ny << ',' << r.seed << ',' << r.flips_requested << ','
          << r.flips_accepted << ',' << r.B << ',' << r.unbounded << ','
          << r.final_E_K << ',' << r.final_max_abs_K << ',' << r.radius_ratio << ','
          << r.polydispersity_delta << ',' << r.boundary_rejections << ',' << r.ricci_iters << ','
          << r.develop_global_vertex_max_spread << ',' << r.period_fit_rms << ','
          << r.nonedge_overlap_count << ',' << r.nonedge_max_overlap << ','
          << r.edge_crossing_count << ',' << r.face_bad_area_count << ','
          << (r.physically_reasonable_radius ? 1 : 0) << ','
          << (r.valid_disk_packing_candidate ? 1 : 0) << ','
          << (r.strong_INC_candidate ? 1 : 0) << ','
          << r.time_total_sec << ',' << '"' << r.status << '"' << '\n';
    }
}

static vector<BoundsRow> run_bounds_scan(Config cfg) {
    apply_output_defaults(cfg, true);
    cfg.progress_every = 0;
    cfg.develop = true;
    if (!cfg.check_crossings_explicit) cfg.check_crossings = false;
    if (cfg.bounds_list.empty()) cfg.bounds_list = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};

    struct Task { double B; int unbounded; };
    vector<Task> tasks;
    if (cfg.include_unbounded) tasks.push_back({numeric_limits<double>::infinity(), 1});
    for (double B : cfg.bounds_list) tasks.push_back({B, 0});

    vector<BoundsRow> rows(tasks.size());
    atomic<size_t> next{0};
    int nthreads = max(1, min(cfg.threads, static_cast<int>(tasks.size())));
    vector<thread> workers;
    for (int t = 0; t < nthreads; ++t) {
        workers.emplace_back([&]() {
            while (true) {
                size_t idx = next.fetch_add(1);
                if (idx >= tasks.size()) break;
                Config local = cfg;
                if (tasks[idx].unbounded) {
                    local.bounded = false;
                } else {
                    local.bounded = true;
                    local.u_min = -tasks[idx].B;
                    local.u_max = tasks[idx].B;
                }
                RunOutput out = execute_run(local, cfg.scan_write_details, false);
                rows[idx] = to_bounds_row(out, tasks[idx].B, tasks[idx].unbounded);
            }
        });
    }
    for (thread& th : workers) th.join();

    for (const BoundsRow& r : rows) {
        cout << scientific << setprecision(6)
             << "bounds B=" << r.B
             << " unbounded=" << r.unbounded
             << " maxK=" << r.final_max_abs_K
             << " ratio=" << r.radius_ratio
             << " overlap=" << r.nonedge_max_overlap
             << " valid=" << (r.valid_disk_packing_candidate ? 1 : 0)
             << " time=" << r.time_total_sec
             << " status=" << r.status << '\n';
    }
    write_bounds_csv(cfg.bounds_out, rows);
    return rows;
}

static double surgery_score(const RunOutput& out, const Config& cfg) {
    double maxK = isfinite(out.ricci.state.max_abs_K) ? out.ricci.state.max_abs_K : 1e100;
    double spread = isfinite(out.develop.global_vertex_max_spread) ? out.develop.global_vertex_max_spread : 1e100;
    double period = isfinite(out.develop.period_fit_rms) ? out.develop.period_fit_rms : 1e100;
    double overlap = isfinite(out.packing.nonedge_max_overlap) ? out.packing.nonedge_max_overlap : 1e100;
    double frac = isfinite(out.packing.nonedge_overlap_fraction) ? out.packing.nonedge_overlap_fraction : 1.0;
    double ratio = isfinite(out.ricci.radius_ratio) && out.ricci.radius_ratio > 0.0 ? out.ricci.radius_ratio : 1e100;
    double delta = isfinite(out.packing.radius_stats.polydispersity_delta) ? out.packing.radius_stats.polydispersity_delta : 1e100;

    if (cfg.surgery_metric == "radius") {
        return 1000.0 * maxK + 1000.0 * spread + 1000.0 * period +
               0.1 * log(max(ratio, 1.0)) + 100.0 * overlap;
    }
    if (cfg.surgery_metric == "mixed") {
        return 1000.0 * maxK + 1000.0 * spread + 1000.0 * period +
               100.0 * overlap + 10.0 * frac +
               0.1 * max(0.0, log(max(ratio / cfg.target_max_ratio, 1e-300))) +
               10.0 * max(0.0, delta - cfg.target_delta);
    }
    return 1000.0 * maxK + 1000.0 * spread + 1000.0 * period +
           100.0 * overlap + 10.0 * frac + 0.01 * log(max(ratio, 1.0));
}

static RunOutput relax_mesh_from_u(Config cfg, Mesh mesh, vector<double> u, bool keep_stats, bool verbose) {
    RunOutput out;
    out.cfg = cfg;
    double t_total0 = now_sec();
    try {
        double t0 = now_sec();
        string error;
        if (!validate_mesh(mesh, &error)) throw runtime_error("mesh invalid: " + error);
        out.mesh = std::move(mesh);
        out.time.topology = now_sec() - t0;

        t0 = now_sec();
        out.ricci = run_ricci(out.mesh, std::move(u), cfg, keep_stats, verbose);
        out.time.ricci = now_sec() - t0;
        if (cfg.develop) {
            t0 = now_sec();
            out.develop = develop_mesh(out.mesh, out.ricci.state.r);
            out.time.develop = now_sec() - t0;
            out.packing = compute_packing_diagnostics(out.mesh, out.ricci, out.develop, cfg);
        }
        out.time.total = now_sec() - t_total0;
    } catch (const exception& e) {
        out.valid = false;
        out.error = e.what();
        out.time.total = now_sec() - t_total0;
    }
    return out;
}

static RunOutput run_surgery(Config cfg) {
    filesystem::create_directories(cfg.out);
    cfg.develop = true;
    cfg.progress_every = 0;
    Config base_cfg = cfg;
    RunOutput current = execute_run(base_cfg, false, false);
    if (!current.valid) return current;
    int initial_flips_accepted = current.flips_accepted;
    int surgery_accepts = 0;
    double current_score = surgery_score(current, cfg);

    ofstream stats(filesystem::path(cfg.out) / "surgery_stats.csv");
    stats << "step,accepted,S_old,S_new,final_max_abs_K,radius_ratio,polydispersity_delta,"
             "develop_spread,period_fit_rms,nonedge_max_overlap,nonedge_overlap_count,edge_crossing_count\n";
    stats << scientific << setprecision(17);

    mt19937_64 rng(cfg.seed + 0x9e3779b97f4a7c15ULL);
    uniform_real_distribution<double> uni(0.0, 1.0);
    for (int step = 1; step <= cfg.surgery_steps; ++step) {
        double old_score = current_score;
        RunOutput best;
        double best_score = numeric_limits<double>::infinity();
        bool have_best = false;

        for (int trial = 0; trial < cfg.surgery_trials_per_step; ++trial) {
            Mesh proposal_mesh = current.mesh;
            if (!attempt_random_flip(proposal_mesh, rng)) continue;
            string error;
            if (!validate_mesh(proposal_mesh, &error)) continue;

            Config trial_cfg = cfg;
            trial_cfg.max_iter = min(cfg.max_iter, cfg.surgery_ricci_iter);
            trial_cfg.progress_every = 0;
            RunOutput proposal = relax_mesh_from_u(trial_cfg, std::move(proposal_mesh), current.ricci.state.u, false, false);
            if (!proposal.valid) continue;
            double S = surgery_score(proposal, cfg);
            if (S < best_score) {
                best_score = S;
                best = std::move(proposal);
                have_best = true;
            }
        }

        bool accepted = false;
        if (have_best && best_score < current_score) {
            accepted = true;
        } else if (have_best && cfg.surgery_temperature > 0.0) {
            double p = exp(-(best_score - current_score) / cfg.surgery_temperature);
            accepted = uni(rng) < p;
        }
        if (accepted) {
            current = std::move(best);
            current_score = best_score;
            surgery_accepts++;
        }

        stats << step << ',' << (accepted ? 1 : 0) << ',' << old_score << ',' << current_score << ','
              << current.ricci.state.max_abs_K << ',' << current.ricci.radius_ratio << ','
              << current.packing.radius_stats.polydispersity_delta << ','
              << current.develop.global_vertex_max_spread << ',' << current.develop.period_fit_rms << ','
              << current.packing.nonedge_max_overlap << ',' << current.packing.nonedge_overlap_count << ','
              << current.packing.edge_crossing_count << '\n';

        if (cfg.surgery_report_every > 0 && step % cfg.surgery_report_every == 0) {
            cout << scientific << setprecision(6)
                 << "surgery step=" << step
                 << " accepted=" << (accepted ? 1 : 0)
                 << " S=" << current_score
                 << " maxK=" << current.ricci.state.max_abs_K
                 << " ratio=" << current.ricci.radius_ratio
                 << " overlap=" << current.packing.nonedge_max_overlap << '\n';
        }
    }
    current.cfg = cfg;
    current.flips_accepted = initial_flips_accepted + surgery_accepts;
    return current;
}

static bool run_tests() {
    bool ok = true;
    cout << "Running v4 tests\n";
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 0; cfg.seed = 1; cfg.max_iter = 1;
        cfg.tol = 1e-12; cfg.develop = true; cfg.progress_every = 0; cfg.check_crossings = true;
        Mesh mesh = generate_initial_mesh(cfg.nx, cfg.ny);
        vector<double> u(mesh.N, 0.0);
        RunOutput out = relax_mesh_from_u(cfg, mesh, u, false, false);
        bool pass = out.valid && out.ricci.state.max_abs_K < 1e-12 &&
                    out.packing.developed_consistent &&
                    out.packing.nonedge_max_overlap <= cfg.overlap_tol &&
                    out.packing.edge_crossing_count == 0;
        ok = ok && pass;
        cout << "Test 1 regular packing diagnostics: " << (pass ? "PASS" : "FAIL")
             << " maxK=" << scientific << setprecision(6) << out.ricci.state.max_abs_K
             << " nonedge_overlap=" << out.packing.nonedge_max_overlap
             << " crossings=" << out.packing.edge_crossing_count
             << " physical=" << (out.packing.physical_INC_like_candidate ? 1 : 0) << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 0; cfg.seed = 2; cfg.max_iter = 100000;
        cfg.tol = 1e-9; cfg.develop = true; cfg.progress_every = 0; cfg.check_crossings = true;
        RunOutput out = execute_run(cfg, false, false);
        bool pass = out.valid && out.ricci.state.max_abs_K < cfg.tol &&
                    out.packing.nonedge_max_overlap <= cfg.overlap_tol &&
                    out.packing.edge_crossing_count == 0;
        ok = ok && pass;
        cout << "Test 2 regular random-u Ricci: " << (pass ? "PASS" : "FAIL")
             << " finalK=" << out.ricci.state.max_abs_K
             << " nonedge_overlap=" << out.packing.nonedge_max_overlap
             << " ratio=" << out.ricci.radius_ratio << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 50; cfg.seed = 1; cfg.max_iter = 100000;
        cfg.tol = 1e-9; cfg.develop = true; cfg.progress_every = 0; cfg.check_crossings = false;
        RunOutput out = execute_run(cfg, false, false);
        bool pass = out.valid && out.develop.ran && out.packing.ran;
        ok = ok && pass;
        cout << "Test 3 random flips unbounded diagnostics: " << (pass ? "PASS" : "FAIL")
             << " finalK=" << out.ricci.state.max_abs_K
             << " nonedge_overlap=" << out.packing.nonedge_max_overlap
             << " valid_disk=" << (out.packing.valid_disk_packing_candidate ? 1 : 0) << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 50; cfg.seed = 1; cfg.max_iter = 100000;
        cfg.tol = 1e-9; cfg.bounded = true; cfg.u_min = -2.0; cfg.u_max = 2.0;
        cfg.develop = true; cfg.progress_every = 0; cfg.check_crossings = false;
        RunOutput out = execute_run(cfg, false, false);
        bool bounds_ok = out.ricci.u_min_actual >= cfg.u_min - 1e-12 && out.ricci.u_max_actual <= cfg.u_max + 1e-12;
        bool pass = out.valid && bounds_ok && out.packing.ran;
        ok = ok && pass;
        cout << "Test 4 random flips bounded diagnostics: " << (pass ? "PASS" : "FAIL")
             << " finalK=" << out.ricci.state.max_abs_K
             << " boundary_rejections=" << out.ricci.boundary_rejections
             << " stopped=" << out.ricci.stopped_reason << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 50; cfg.seed = 1; cfg.max_iter = 50000;
        cfg.tol = 1e-9; cfg.develop = true; cfg.scan_bounds = true; cfg.bounds_list = {0.5, 1.0, 2.0, 3.0};
        cfg.bounds_out = "bounds_scan_test.csv"; cfg.progress_every = 0; cfg.scan_write_details = false;
        vector<BoundsRow> rows = run_bounds_scan(cfg);
        bool file_exists = filesystem::exists(cfg.bounds_out);
        bool pass = rows.size() == 4 && file_exists;
        ok = ok && pass;
        cout << "Test 5 bounds scan small: " << (pass ? "PASS" : "FAIL")
             << " rows=" << rows.size()
             << " file=" << (file_exists ? 1 : 0) << '\n';
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
        if (cfg.scan_bounds) {
            double t0 = now_sec();
            vector<BoundsRow> rows = run_bounds_scan(cfg);
            double total = now_sec() - t0;
            cout << "bounds_rows " << rows.size() << " bounds_out " << cfg.bounds_out
                 << " total_sec " << scientific << setprecision(6) << total
                 << " threads " << cfg.threads << '\n';
            return 0;
        }

        RunOutput out = cfg.surgery ? run_surgery(cfg) : execute_run(cfg, true, cfg.progress_every > 0);
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
        cout << "polydispersity_delta " << out.packing.radius_stats.polydispersity_delta << '\n';
        cout << "max_edge_contact_error " << out.packing.max_edge_contact_error << '\n';
        cout << "nonedge_overlap_count " << out.packing.nonedge_overlap_count << '\n';
        cout << "nonedge_max_overlap " << out.packing.nonedge_max_overlap << '\n';
        cout << "edge_crossing_count " << out.packing.edge_crossing_count << '\n';
        cout << "face_bad_area_count " << out.packing.face_bad_area_count << '\n';
        cout << "valid_disk_packing_candidate " << (out.packing.valid_disk_packing_candidate ? "yes" : "no") << '\n';
        cout << "strong_INC_candidate " << (out.packing.strong_INC_candidate ? "yes" : "no") << '\n';
        cout << "physical_INC_like_candidate " << (out.packing.physical_INC_like_candidate ? "yes" : "no") << '\n';
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
