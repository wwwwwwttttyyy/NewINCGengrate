#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
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

struct EdgeKey {
    int u, v;

    EdgeKey() : u(0), v(0) {}
    EdgeKey(int x, int y) {
        if (x < y) {
            u = x;
            v = y;
        } else {
            u = y;
            v = x;
        }
    }

    bool operator<(const EdgeKey& other) const {
        if (u != other.u) return u < other.u;
        return v < other.v;
    }
};

struct EdgeUse {
    int face = -1;
    int p = -1;
    int q = -1;
};

struct EdgeInfo {
    int i = 0;
    int j = 0;
    int shift_x = 0;
    int shift_y = 0;
};

struct TriangulationStats {
    int V = 0;
    int E = 0;
    int F = 0;
    int chi = 0;
};

struct State {
    vector<double> r;
    vector<double> K;
    vector<array<double, 3>> face_angles;
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

struct EmbedStatsRow {
    int iter = 0;
    double E_embed = 0.0;
    double max_edge_abs_residual = 0.0;
    double rms_edge_residual = 0.0;
    double Lx = 0.0;
    double Ly = 0.0;
    double shear = 0.0;
};

struct RunResult {
    vector<double> u;
    State state;
    vector<StatRow> stats;
    int iterations = 0;
    string stopped_reason;
    bool monotonic_E = true;
    double initial_E_K = 0.0;
    double initial_max_abs_K = 0.0;
    int boundary_rejections = 0;
    double u_min_actual = 0.0;
    double u_max_actual = 0.0;
    double radius_ratio = 0.0;
};

struct EmbedResult {
    bool ran = false;
    vector<double> sx;
    vector<double> sy;
    vector<EmbedStatsRow> stats;
    int iterations = 0;
    string stopped_reason = "not requested";
    double E_embed = numeric_limits<double>::quiet_NaN();
    double max_edge_abs_residual = numeric_limits<double>::quiet_NaN();
    double rms_edge_residual = numeric_limits<double>::quiet_NaN();
    double Lx = numeric_limits<double>::quiet_NaN();
    double Ly = numeric_limits<double>::quiet_NaN();
    double shear = numeric_limits<double>::quiet_NaN();
};

struct ThetaResult {
    vector<double> theta_vertex;
    vector<int> theta_vertex_count;
    double theta_mean_abs = numeric_limits<double>::quiet_NaN();
    double theta_rms = numeric_limits<double>::quiet_NaN();
    double theta_max_abs = numeric_limits<double>::quiet_NaN();
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

    bool embed = false;
    int embed_iter = 200000;
    double embed_tol = 1e-10;
    double embed_lr = 0.01;
};

struct LocalFaceRef {
    int face = -1;
    array<int, 3> affected_pos{};
};

struct LocalStencil {
    vector<int> affected_vertices;
    vector<LocalFaceRef> incident_faces;
};

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

static int floor_div_exact_periodic(int x, int n) {
    int base = imod(x, n);
    return (x - base) / n;
}

static double wrap01(double x) {
    x -= floor(x);
    if (x >= 1.0) x -= 1.0;
    if (x < 0.0) x += 1.0;
    return x;
}

static Corner corner_from_unwrapped(int i, int j, int nx, int ny) {
    return {vertex_id(i, j, nx, ny), floor_div_exact_periodic(i, nx), floor_div_exact_periodic(j, ny)};
}

static Face make_face(Corner a, Corner b, Corner c) {
    Face f;
    f.c[0] = a;
    f.c[1] = b;
    f.c[2] = c;
    return f;
}

static vector<Face> generate_periodic_triangular_torus(int nx, int ny) {
    vector<Face> faces;
    faces.reserve(static_cast<size_t>(2 * nx * ny));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            Corner A = corner_from_unwrapped(i, j, nx, ny);
            Corner B = corner_from_unwrapped(i + 1, j, nx, ny);
            Corner C = corner_from_unwrapped(i, j + 1, nx, ny);
            Corner D = corner_from_unwrapped(i + 1, j + 1, nx, ny);
            faces.push_back(make_face(A, B, C));
            faces.push_back(make_face(B, D, C));
        }
    }
    return faces;
}

static map<EdgeKey, vector<EdgeUse>> build_edges(const vector<Face>& faces) {
    map<EdgeKey, vector<EdgeUse>> edge_to_faces;
    for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi) {
        const Face& f = faces[fi];
        const int p[3] = {0, 1, 2};
        const int q[3] = {1, 2, 0};
        for (int e = 0; e < 3; ++e) {
            edge_to_faces[EdgeKey(f.c[p[e]].v, f.c[q[e]].v)].push_back({fi, p[e], q[e]});
        }
    }
    return edge_to_faces;
}

static bool canonical_shift_for_use(const Face& f, const EdgeUse& use, int* sx, int* sy) {
    const Corner& a = f.c[use.p];
    const Corner& b = f.c[use.q];
    EdgeKey key(a.v, b.v);
    int dx = b.sx - a.sx;
    int dy = b.sy - a.sy;
    if (a.v == key.u && b.v == key.v) {
        *sx = dx;
        *sy = dy;
        return true;
    }
    if (a.v == key.v && b.v == key.u) {
        *sx = -dx;
        *sy = -dy;
        return true;
    }
    return false;
}

static vector<EdgeInfo> build_edge_infos(const vector<Face>& faces, string* error = nullptr) {
    auto edge_to_faces = build_edges(faces);
    vector<EdgeInfo> infos;
    infos.reserve(edge_to_faces.size());

    for (const auto& kv : edge_to_faces) {
        bool have_shift = false;
        int sx0 = 0;
        int sy0 = 0;
        for (const EdgeUse& use : kv.second) {
            int sx = 0;
            int sy = 0;
            if (!canonical_shift_for_use(faces[use.face], use, &sx, &sy)) {
                if (error) *error = "failed to compute canonical edge shift";
                return {};
            }
            if (!have_shift) {
                sx0 = sx;
                sy0 = sy;
                have_shift = true;
            } else if (sx != sx0 || sy != sy0) {
                if (error) {
                    ostringstream oss;
                    oss << "inconsistent periodic shift on edge (" << kv.first.u << "," << kv.first.v << ")";
                    *error = oss.str();
                }
                return {};
            }
        }
        infos.push_back({kv.first.u, kv.first.v, sx0, sy0});
    }
    return infos;
}

static bool validate_triangulation(const vector<Face>& faces, int V, TriangulationStats* out_stats = nullptr,
                                   string* error = nullptr) {
    auto fail = [&](const string& msg) {
        if (error) *error = msg;
        return false;
    };

    for (const Face& f : faces) {
        for (int q = 0; q < 3; ++q) {
            if (f.c[q].v < 0 || f.c[q].v >= V) return fail("face has out-of-range vertex index");
        }
        if (f.c[0].v == f.c[1].v || f.c[1].v == f.c[2].v || f.c[2].v == f.c[0].v) {
            return fail("face has repeated vertices");
        }
    }

    auto edge_to_faces = build_edges(faces);
    for (const auto& kv : edge_to_faces) {
        if (kv.second.size() != 2) return fail("edge does not have exactly two incident faces");
    }

    string shift_error;
    build_edge_infos(faces, &shift_error);
    if (!shift_error.empty()) return fail(shift_error);

    TriangulationStats stats;
    stats.V = V;
    stats.E = static_cast<int>(edge_to_faces.size());
    stats.F = static_cast<int>(faces.size());
    stats.chi = stats.V - stats.E + stats.F;
    if (out_stats) *out_stats = stats;

    if (stats.chi != 0) return fail("Euler characteristic is not zero");
    if (stats.F != 2 * stats.V) return fail("face count is not 2V");
    if (stats.E != 3 * stats.V) return fail("edge count is not 3V");
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

static bool nondegenerate(const Face& f) {
    return f.c[0].v != f.c[1].v && f.c[1].v != f.c[2].v && f.c[2].v != f.c[0].v;
}

static Corner shifted_corner(Corner c, int ox, int oy) {
    c.sx += ox;
    c.sy += oy;
    return c;
}

static bool attempt_random_flip(vector<Face>& faces, int V, mt19937_64& rng) {
    auto edge_to_faces = build_edges(faces);
    if (edge_to_faces.empty()) return false;

    vector<EdgeKey> edges;
    edges.reserve(edge_to_faces.size());
    for (const auto& kv : edge_to_faces) edges.push_back(kv.first);

    uniform_int_distribution<size_t> pick(0, edges.size() - 1);
    const EdgeKey e = edges[pick(rng)];
    const auto& inc = edge_to_faces[e];
    if (inc.size() != 2) return false;

    int f1 = inc[0].face;
    int f2 = inc[1].face;
    const Face& face1 = faces[f1];
    const Face& face2 = faces[f2];
    int iu1 = corner_index_for_vertex(face1, e.u);
    int iv1 = corner_index_for_vertex(face1, e.v);
    int iu2 = corner_index_for_vertex(face2, e.u);
    int iv2 = corner_index_for_vertex(face2, e.v);
    int io1 = opposite_index(face1, e.u, e.v);
    int io2 = opposite_index(face2, e.u, e.v);
    if (iu1 < 0 || iv1 < 0 || iu2 < 0 || iv2 < 0 || io1 < 0 || io2 < 0) return false;

    Corner a = face1.c[iu1];
    Corner b = face1.c[iv1];
    Corner c = face1.c[io1];

    int ox = a.sx - face2.c[iu2].sx;
    int oy = a.sy - face2.c[iu2].sy;
    Corner b2 = shifted_corner(face2.c[iv2], ox, oy);
    if (b2.sx != b.sx || b2.sy != b.sy) return false;
    Corner d = shifted_corner(face2.c[io2], ox, oy);

    if (c.v == d.v) return false;
    if (edge_to_faces.find(EdgeKey(c.v, d.v)) != edge_to_faces.end()) return false;

    Face nf1 = make_face(c, d, a);
    Face nf2 = make_face(d, c, b);
    if (!nondegenerate(nf1) || !nondegenerate(nf2)) return false;

    Face old1 = faces[f1];
    Face old2 = faces[f2];
    faces[f1] = nf1;
    faces[f2] = nf2;

    if (!validate_triangulation(faces, V)) {
        faces[f1] = old1;
        faces[f2] = old2;
        return false;
    }
    return true;
}

static void normalize_u(vector<double>& u) {
    if (u.empty()) return;
    double mean = accumulate(u.begin(), u.end(), 0.0) / static_cast<double>(u.size());
    for (double& x : u) x -= mean;
}

static vector<double> random_log_radii(int N, mt19937_64& rng, double stddev = 0.1) {
    normal_distribution<double> normal(0.0, stddev);
    vector<double> u(N);
    for (double& x : u) x = normal(rng);
    normalize_u(u);
    return u;
}

static double clamp_unit(double x) {
    if (x < -1.0) return -1.0;
    if (x > 1.0) return 1.0;
    return x;
}

static array<double, 3> compute_angles_from_radii(double ra, double rb, double rc) {
    auto angle_at = [](double ri, double rj, double rk) {
        double x = ri + rj;
        double y = ri + rk;
        double z = rj + rk;
        double scale = max(x, max(y, z));
        if (scale <= 0.0 || !isfinite(scale)) return numeric_limits<double>::quiet_NaN();
        x /= scale;
        y /= scale;
        z /= scale;
        double denom = 2.0 * x * y;
        if (denom <= 0.0) return numeric_limits<double>::quiet_NaN();
        return acos(clamp_unit((x * x + y * y - z * z) / denom));
    };

    return {angle_at(ra, rb, rc), angle_at(rb, rc, ra), angle_at(rc, ra, rb)};
}

static array<double, 3> compute_angles_for_face(const Face& f, const vector<double>& r) {
    return compute_angles_from_radii(r[f.c[0].v], r[f.c[1].v], r[f.c[2].v]);
}

static double safe_exp(double x) {
    if (x > 700.0) x = 700.0;
    if (x < -700.0) x = -700.0;
    return exp(x);
}

static State compute_state(const vector<double>& u, const vector<Face>& faces) {
    State s;
    const int N = static_cast<int>(u.size());
    s.r.resize(N);
    s.K.assign(N, 2.0 * PI);
    s.face_angles.resize(faces.size());
    s.r_min = numeric_limits<double>::infinity();
    s.r_max = 0.0;

    for (int i = 0; i < N; ++i) {
        s.r[i] = safe_exp(u[i]);
        s.r_min = min(s.r_min, s.r[i]);
        s.r_max = max(s.r_max, s.r[i]);
    }

    for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi) {
        const Face& f = faces[fi];
        auto ang = compute_angles_for_face(f, s.r);
        s.face_angles[fi] = ang;
        s.K[f.c[0].v] -= ang[0];
        s.K[f.c[1].v] -= ang[1];
        s.K[f.c[2].v] -= ang[2];
    }

    s.E_K = 0.0;
    s.max_abs_K = 0.0;
    for (double k : s.K) {
        s.E_K += 0.5 * k * k;
        s.max_abs_K = max(s.max_abs_K, fabs(k));
    }
    return s;
}

static vector<double> compute_curvature(const vector<double>& u, const vector<Face>& faces) {
    return compute_state(u, faces).K;
}

static StatRow make_stat_row(int iter, const State& s, double step, int accepted) {
    return {iter, s.E_K, s.max_abs_K, step, s.r_min, s.r_max, accepted};
}

static void print_progress(const StatRow& row) {
    cout << scientific << setprecision(6)
         << "iter " << row.iter
         << " E_K " << row.E_K
         << " max_abs_K " << row.max_abs_K
         << " step " << row.step
         << " r_min " << row.r_min
         << " r_max " << row.r_max
         << " accepted " << row.accepted << '\n';
}

static void finalize_run_diagnostics(RunResult& result) {
    if (!result.u.empty()) {
        auto mm = minmax_element(result.u.begin(), result.u.end());
        result.u_min_actual = *mm.first;
        result.u_max_actual = *mm.second;
    }
    result.radius_ratio = result.state.r_min > 0.0 ? result.state.r_max / result.state.r_min
                                                   : numeric_limits<double>::infinity();
}

static bool inside_u_bounds(const vector<double>& u, double u_min, double u_max) {
    for (double x : u) {
        if (x < u_min || x > u_max) return false;
    }
    return true;
}

static vector<LocalStencil> build_local_stencils(const vector<Face>& faces, int N) {
    vector<vector<int>> incident(N);
    for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi) {
        const Face& f = faces[fi];
        for (int q = 0; q < 3; ++q) incident[f.c[q].v].push_back(fi);
    }

    vector<LocalStencil> stencils(N);
    for (int i = 0; i < N; ++i) {
        vector<int> affected;
        for (int fi : incident[i]) {
            const Face& f = faces[fi];
            for (int q = 0; q < 3; ++q) {
                int v = f.c[q].v;
                if (find(affected.begin(), affected.end(), v) == affected.end()) affected.push_back(v);
            }
        }
        stencils[i].affected_vertices = affected;
        for (int fi : incident[i]) {
            const Face& f = faces[fi];
            LocalFaceRef ref;
            ref.face = fi;
            for (int q = 0; q < 3; ++q) {
                auto it = find(affected.begin(), affected.end(), f.c[q].v);
                ref.affected_pos[q] = static_cast<int>(it - affected.begin());
            }
            stencils[i].incident_faces.push_back(ref);
        }
    }
    return stencils;
}

static vector<double> finite_difference_gradient(const vector<double>& u, const vector<Face>& faces,
                                                 const vector<LocalStencil>& stencils, const State& state,
                                                 double h) {
    const int N = static_cast<int>(u.size());
    vector<double> grad(N, 0.0);
    const double eh = exp(h);
    const double emh = exp(-h);

    for (int i = 0; i < N; ++i) {
        const LocalStencil& stencil = stencils[i];
        const int M = static_cast<int>(stencil.affected_vertices.size());
        vector<double> k_plus(M), k_minus(M);
        double old_local_E = 0.0;
        for (int p = 0; p < M; ++p) {
            int v = stencil.affected_vertices[p];
            k_plus[p] = state.K[v];
            k_minus[p] = state.K[v];
            old_local_E += 0.5 * state.K[v] * state.K[v];
        }

        const double ri_plus = state.r[i] * eh;
        const double ri_minus = state.r[i] * emh;

        for (const LocalFaceRef& ref : stencil.incident_faces) {
            const Face& f = faces[ref.face];
            double rp[3] = {state.r[f.c[0].v], state.r[f.c[1].v], state.r[f.c[2].v]};
            double rm[3] = {state.r[f.c[0].v], state.r[f.c[1].v], state.r[f.c[2].v]};
            for (int q = 0; q < 3; ++q) {
                if (f.c[q].v == i) {
                    rp[q] = ri_plus;
                    rm[q] = ri_minus;
                }
            }

            auto plus_angles = compute_angles_from_radii(rp[0], rp[1], rp[2]);
            auto minus_angles = compute_angles_from_radii(rm[0], rm[1], rm[2]);
            const auto& old_angles = state.face_angles[ref.face];
            for (int q = 0; q < 3; ++q) {
                int p = ref.affected_pos[q];
                k_plus[p] += old_angles[q] - plus_angles[q];
                k_minus[p] += old_angles[q] - minus_angles[q];
            }
        }

        double e_plus = state.E_K - old_local_E;
        double e_minus = state.E_K - old_local_E;
        for (int p = 0; p < M; ++p) {
            e_plus += 0.5 * k_plus[p] * k_plus[p];
            e_minus += 0.5 * k_minus[p] * k_minus[p];
        }
        grad[i] = (e_plus - e_minus) / (2.0 * h);
    }

    double mean = accumulate(grad.begin(), grad.end(), 0.0) / static_cast<double>(N);
    for (double& g : grad) g -= mean;
    return grad;
}

static RunResult run_ricci(vector<double> u, const vector<Face>& faces, int max_iter, double tol,
                           bool bounded, double u_min, double u_max, bool verbose) {
    normalize_u(u);
    RunResult result;
    result.u = u;

    double dt = 0.01;
    const double min_step = 1e-14;
    State state = compute_state(result.u, faces);
    result.initial_E_K = state.E_K;
    result.initial_max_abs_K = state.max_abs_K;
    result.stats.push_back(make_stat_row(0, state, dt, 1));
    if (verbose) print_progress(result.stats.back());

    if (state.max_abs_K < tol) {
        result.state = state;
        result.iterations = 0;
        result.stopped_reason = "converged: max_abs_K < tol at iter 0";
        finalize_run_diagnostics(result);
        return result;
    }

    for (int iter = 1; iter <= max_iter; ++iter) {
        bool accepted = false;
        double accepted_step = dt;
        State candidate_state;
        vector<double> candidate_u(result.u.size());

        while (!accepted) {
            if (dt < min_step) {
                result.state = state;
                result.iterations = iter - 1;
                result.stopped_reason = "stalled: dt < 1e-14";
                result.stats.push_back(make_stat_row(result.iterations, state, dt, 0));
                finalize_run_diagnostics(result);
                return result;
            }

            for (size_t i = 0; i < result.u.size(); ++i) candidate_u[i] = result.u[i] - dt * state.K[i];
            normalize_u(candidate_u);

            if (bounded && !inside_u_bounds(candidate_u, u_min, u_max)) {
                result.boundary_rejections++;
                dt *= 0.5;
                continue;
            }

            candidate_state = compute_state(candidate_u, faces);
            if (isfinite(candidate_state.E_K) && candidate_state.E_K < state.E_K) {
                accepted = true;
                accepted_step = dt;
                result.monotonic_E = result.monotonic_E && (candidate_state.E_K <= state.E_K + 1e-15);
                result.u.swap(candidate_u);
                state = std::move(candidate_state);
                dt = min(dt * 1.05, 1.0);
            } else {
                dt *= 0.5;
            }
        }

        result.iterations = iter;
        if (iter % 100 == 0 || state.max_abs_K < tol) {
            result.stats.push_back(make_stat_row(iter, state, accepted_step, 1));
            if (verbose) print_progress(result.stats.back());
        }

        if (state.max_abs_K < tol) {
            result.state = state;
            result.stopped_reason = "converged: max_abs_K < tol";
            finalize_run_diagnostics(result);
            return result;
        }
    }

    result.state = state;
    result.stopped_reason = "reached max_iter";
    if (result.stats.empty() || result.stats.back().iter != result.iterations) {
        result.stats.push_back(make_stat_row(result.iterations, result.state, dt, 1));
    }
    finalize_run_diagnostics(result);
    return result;
}

static RunResult run_grad(vector<double> u, const vector<Face>& faces, int max_iter, double tol, bool verbose) {
    normalize_u(u);
    RunResult result;
    result.u = u;

    const int N = static_cast<int>(u.size());
    const double h = 1e-6;
    const double alpha_start = 0.1;
    double last_alpha = alpha_start;
    const double armijo_c = 1e-4;
    const double min_step = 1e-14;
    vector<LocalStencil> stencils = build_local_stencils(faces, N);

    State state = compute_state(result.u, faces);
    result.initial_E_K = state.E_K;
    result.initial_max_abs_K = state.max_abs_K;
    result.stats.push_back(make_stat_row(0, state, last_alpha, 1));
    if (verbose) print_progress(result.stats.back());

    if (state.max_abs_K < tol) {
        result.state = state;
        result.iterations = 0;
        result.stopped_reason = "converged: max_abs_K < tol at iter 0";
        finalize_run_diagnostics(result);
        return result;
    }

    for (int iter = 1; iter <= max_iter; ++iter) {
        vector<double> grad = finite_difference_gradient(result.u, faces, stencils, state, h);
        double grad_norm_sq = 0.0;
        for (double g : grad) grad_norm_sq += g * g;
        if (grad_norm_sq < 1e-28) {
            result.state = state;
            result.iterations = iter - 1;
            result.stopped_reason = "stalled: projected gradient norm too small";
            result.stats.push_back(make_stat_row(result.iterations, state, 0.0, 0));
            finalize_run_diagnostics(result);
            return result;
        }

        bool accepted = false;
        double trial_alpha = alpha_start;
        double accepted_step = trial_alpha;
        State candidate_state;
        vector<double> candidate_u(result.u.size());

        while (!accepted) {
            if (trial_alpha < min_step) {
                result.state = state;
                result.iterations = iter - 1;
                result.stopped_reason = "stalled: alpha < 1e-14";
                result.stats.push_back(make_stat_row(result.iterations, state, trial_alpha, 0));
                finalize_run_diagnostics(result);
                return result;
            }

            for (size_t i = 0; i < result.u.size(); ++i) candidate_u[i] = result.u[i] - trial_alpha * grad[i];
            normalize_u(candidate_u);
            candidate_state = compute_state(candidate_u, faces);

            double armijo_rhs = state.E_K - armijo_c * trial_alpha * grad_norm_sq;
            if (isfinite(candidate_state.E_K) && candidate_state.E_K <= armijo_rhs) {
                accepted = true;
                accepted_step = trial_alpha;
                result.monotonic_E = result.monotonic_E && (candidate_state.E_K <= state.E_K + 1e-15);
                result.u.swap(candidate_u);
                state = std::move(candidate_state);
                last_alpha = accepted_step;
            } else {
                trial_alpha *= 0.5;
            }
        }

        result.iterations = iter;
        if (iter % 100 == 0 || state.max_abs_K < tol) {
            result.stats.push_back(make_stat_row(iter, state, accepted_step, 1));
            if (verbose) print_progress(result.stats.back());
        }

        if (state.max_abs_K < tol) {
            result.state = state;
            result.stopped_reason = "converged: max_abs_K < tol";
            finalize_run_diagnostics(result);
            return result;
        }
    }

    result.state = state;
    result.stopped_reason = "reached max_iter";
    if (result.stats.empty() || result.stats.back().iter != result.iterations) {
        result.stats.push_back(make_stat_row(result.iterations, result.state, last_alpha, 1));
    }
    finalize_run_diagnostics(result);
    return result;
}

struct EmbedEval {
    double E = 0.0;
    double max_abs = 0.0;
    double rms = 0.0;
    vector<double> gx;
    vector<double> gy;
    double gLx = 0.0;
    double gLy = 0.0;
    double gShear = 0.0;
};

static EmbedEval evaluate_embedding(const vector<EdgeInfo>& edges, const vector<double>& r,
                                    const vector<double>& sx, const vector<double>& sy,
                                    double Lx, double Ly, double shear, bool gradients) {
    const int N = static_cast<int>(sx.size());
    EmbedEval ev;
    if (gradients) {
        ev.gx.assign(N, 0.0);
        ev.gy.assign(N, 0.0);
    }

    for (const EdgeInfo& e : edges) {
        double dx = sx[e.j] + static_cast<double>(e.shift_x) - sx[e.i];
        double dy = sy[e.j] + static_cast<double>(e.shift_y) - sy[e.i];
        double vx = Lx * dx + shear * dy;
        double vy = Ly * dy;
        double dist = sqrt(vx * vx + vy * vy);
        double target = r[e.i] + r[e.j];
        double residual = dist - target;
        ev.E += 0.5 * residual * residual;
        ev.max_abs = max(ev.max_abs, fabs(residual));
        ev.rms += residual * residual;

        if (gradients) {
            double inv = 1.0 / max(dist, 1e-12);
            double gvx = residual * vx * inv;
            double gvy = residual * vy * inv;
            double gdx = Lx * gvx;
            double gdy = shear * gvx + Ly * gvy;
            ev.gx[e.i] -= gdx;
            ev.gy[e.i] -= gdy;
            ev.gx[e.j] += gdx;
            ev.gy[e.j] += gdy;
            ev.gLx += gvx * dx;
            ev.gShear += gvx * dy;
            ev.gLy += gvy * dy;
        }
    }

    if (!edges.empty()) ev.rms = sqrt(ev.rms / static_cast<double>(edges.size()));
    return ev;
}

static EmbedStatsRow make_embed_row(int iter, const EmbedEval& ev, double Lx, double Ly, double shear) {
    return {iter, ev.E, ev.max_abs, ev.rms, Lx, Ly, shear};
}

static void initialize_fractional_coords(int nx, int ny, vector<double>& sx, vector<double>& sy) {
    const int N = nx * ny;
    sx.assign(N, 0.0);
    sy.assign(N, 0.0);
    for (int v = 0; v < N; ++v) {
        sx[v] = static_cast<double>(vertex_ix(v, nx)) / static_cast<double>(nx);
        sy[v] = static_cast<double>(vertex_iy(v, nx)) / static_cast<double>(ny);
    }
}

static bool solve_3x3(double A[3][3], double b[3], double x[3]) {
    double aug[3][4];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) aug[i][j] = A[i][j];
        aug[i][3] = b[i];
    }

    for (int col = 0; col < 3; ++col) {
        int piv = col;
        for (int r = col + 1; r < 3; ++r) {
            if (fabs(aug[r][col]) > fabs(aug[piv][col])) piv = r;
        }
        if (fabs(aug[piv][col]) < 1e-18) return false;
        if (piv != col) {
            for (int c = col; c < 4; ++c) swap(aug[piv][c], aug[col][c]);
        }
        double div = aug[col][col];
        for (int c = col; c < 4; ++c) aug[col][c] /= div;
        for (int r = 0; r < 3; ++r) {
            if (r == col) continue;
            double factor = aug[r][col];
            for (int c = col; c < 4; ++c) aug[r][c] -= factor * aug[col][c];
        }
    }

    for (int i = 0; i < 3; ++i) x[i] = aug[i][3];
    return true;
}

static bool metric_warm_start(const vector<EdgeInfo>& edges, const vector<double>& r,
                              const vector<double>& sx, const vector<double>& sy,
                              double& Lx, double& Ly, double& shear) {
    double M[3][3] = {};
    double rhs[3] = {};
    for (const EdgeInfo& e : edges) {
        double dx = sx[e.j] + static_cast<double>(e.shift_x) - sx[e.i];
        double dy = sy[e.j] + static_cast<double>(e.shift_y) - sy[e.i];
        double coeff[3] = {dx * dx, 2.0 * dx * dy, dy * dy};
        double target2 = (r[e.i] + r[e.j]) * (r[e.i] + r[e.j]);
        for (int a = 0; a < 3; ++a) {
            rhs[a] += coeff[a] * target2;
            for (int b = 0; b < 3; ++b) M[a][b] += coeff[a] * coeff[b];
        }
    }

    double g[3];
    if (!solve_3x3(M, rhs, g)) return false;
    double g11 = g[0];
    double g12 = g[1];
    double g22 = g[2];
    if (!(g11 > 1e-18)) return false;
    double cand_Lx = sqrt(g11);
    double cand_shear = g12 / cand_Lx;
    double ly2 = g22 - cand_shear * cand_shear;
    if (!(ly2 > 1e-18)) return false;
    double cand_Ly = sqrt(ly2);
    if (!isfinite(cand_Lx) || !isfinite(cand_Ly) || !isfinite(cand_shear)) return false;

    EmbedEval old_ev = evaluate_embedding(edges, r, sx, sy, Lx, Ly, shear, false);
    EmbedEval new_ev = evaluate_embedding(edges, r, sx, sy, cand_Lx, cand_Ly, cand_shear, false);
    if (new_ev.E <= old_ev.E) {
        Lx = cand_Lx;
        Ly = cand_Ly;
        shear = cand_shear;
        return true;
    }
    return false;
}

static EmbedResult run_embedding(const vector<EdgeInfo>& edges, const vector<double>& r, int nx, int ny,
                                 int max_iter, double tol, double lr, bool verbose = false) {
    EmbedResult result;
    result.ran = true;
    initialize_fractional_coords(nx, ny, result.sx, result.sy);

    result.Lx = static_cast<double>(nx);
    result.Ly = static_cast<double>(ny);
    result.shear = 0.0;
    metric_warm_start(edges, r, result.sx, result.sy, result.Lx, result.Ly, result.shear);
    const double fixed_x0 = result.sx[0];
    const double fixed_y0 = result.sy[0];

    const int N = nx * ny;
    vector<double> mx(N, 0.0), my(N, 0.0), vx(N, 0.0), vy(N, 0.0);
    double mLx = 0.0, mLy = 0.0, mSh = 0.0;
    double vLx = 0.0, vLy = 0.0, vSh = 0.0;
    const double beta1 = 0.9;
    const double beta2 = 0.999;
    const double eps = 1e-8;

    auto adam_update = [&](double grad, double& value, double& m, double& v, int iter, double step_scale = 1.0) {
        m = beta1 * m + (1.0 - beta1) * grad;
        v = beta2 * v + (1.0 - beta2) * grad * grad;
        double mhat = m / (1.0 - pow(beta1, iter));
        double vhat = v / (1.0 - pow(beta2, iter));
        value -= lr * step_scale * mhat / (sqrt(vhat) + eps);
    };

    EmbedEval initial = evaluate_embedding(edges, r, result.sx, result.sy, result.Lx, result.Ly, result.shear, false);
    result.stats.push_back(make_embed_row(0, initial, result.Lx, result.Ly, result.shear));
    if (initial.max_abs < tol) {
        result.E_embed = initial.E;
        result.max_edge_abs_residual = initial.max_abs;
        result.rms_edge_residual = initial.rms;
        result.iterations = 0;
        result.stopped_reason = "converged: max_edge_abs_residual < embed_tol at iter 0";
        return result;
    }

    EmbedEval current = initial;
    for (int iter = 1; iter <= max_iter; ++iter) {
        EmbedEval ev = evaluate_embedding(edges, r, result.sx, result.sy, result.Lx, result.Ly, result.shear, true);

        for (int i = 1; i < N; ++i) {
            adam_update(ev.gx[i], result.sx[i], mx[i], vx[i], iter);
            adam_update(ev.gy[i], result.sy[i], my[i], vy[i], iter);
        }

        adam_update(ev.gLx, result.Lx, mLx, vLx, iter, 10.0);
        adam_update(ev.gLy, result.Ly, mLy, vLy, iter, 10.0);
        adam_update(ev.gShear, result.shear, mSh, vSh, iter, 10.0);

        result.Lx = max(result.Lx, 1e-8);
        result.Ly = max(result.Ly, 1e-8);
        result.sx[0] = fixed_x0;
        result.sy[0] = fixed_y0;
        for (int i = 1; i < N; ++i) {
            result.sx[i] = wrap01(result.sx[i]);
            result.sy[i] = wrap01(result.sy[i]);
        }

        current = evaluate_embedding(edges, r, result.sx, result.sy, result.Lx, result.Ly, result.shear, false);
        result.iterations = iter;
        if (iter % 100 == 0 || current.max_abs < tol) {
            result.stats.push_back(make_embed_row(iter, current, result.Lx, result.Ly, result.shear));
        }
        if (verbose && (iter % 10000 == 0 || current.max_abs < tol)) {
            cout << scientific << setprecision(6)
                 << "embed_iter " << iter
                 << " E_embed " << current.E
                 << " max_edge_abs_residual " << current.max_abs
                 << " rms_edge_residual " << current.rms
                 << " Lx " << result.Lx
                 << " Ly " << result.Ly
                 << " shear " << result.shear << '\n';
        }
        if (current.max_abs < tol) {
            result.E_embed = current.E;
            result.max_edge_abs_residual = current.max_abs;
            result.rms_edge_residual = current.rms;
            result.stopped_reason = "converged: max_edge_abs_residual < embed_tol";
            return result;
        }
    }

    result.E_embed = current.E;
    result.max_edge_abs_residual = current.max_abs;
    result.rms_edge_residual = current.rms;
    result.stopped_reason = "reached embed_iter";
    if (result.stats.empty() || result.stats.back().iter != result.iterations) {
        result.stats.push_back(make_embed_row(result.iterations, current, result.Lx, result.Ly, result.shear));
    }
    return result;
}

static double angle_between(double ax, double ay, double bx, double by) {
    double na = sqrt(ax * ax + ay * ay);
    double nb = sqrt(bx * bx + by * by);
    if (na <= 0.0 || nb <= 0.0) return numeric_limits<double>::quiet_NaN();
    return acos(clamp_unit((ax * bx + ay * by) / (na * nb)));
}

static ThetaResult compute_theta_residual(const vector<Face>& faces, const vector<double>& r,
                                          const EmbedResult& embed) {
    const int N = static_cast<int>(r.size());
    ThetaResult tr;
    tr.theta_vertex.assign(N, 0.0);
    tr.theta_vertex_count.assign(N, 0);

    double sum_abs = 0.0;
    double sum_sq = 0.0;
    double max_abs = 0.0;
    int count = 0;

    for (const Face& f : faces) {
        double px[3], py[3];
        for (int q = 0; q < 3; ++q) {
            double fx = embed.sx[f.c[q].v] + static_cast<double>(f.c[q].sx);
            double fy = embed.sy[f.c[q].v] + static_cast<double>(f.c[q].sy);
            px[q] = embed.Lx * fx + embed.shear * fy;
            py[q] = embed.Ly * fy;
        }

        for (int q = 0; q < 3; ++q) {
            int q1 = (q + 1) % 3;
            int q2 = (q + 2) % 3;
            double actual = angle_between(px[q1] - px[q], py[q1] - py[q],
                                          px[q2] - px[q], py[q2] - py[q]);
            double ri = r[f.c[q].v];
            double rj = r[f.c[q1].v];
            double rk = r[f.c[q2].v];
            double a = ri + rj;
            double b = ri + rk;
            double c = rj + rk;
            double ideal = acos(clamp_unit((a * a + b * b - c * c) / (2.0 * a * b)));
            double err = fabs(actual - ideal);
            if (!isfinite(err)) continue;
            tr.theta_vertex[f.c[q].v] += err;
            tr.theta_vertex_count[f.c[q].v]++;
            sum_abs += err;
            sum_sq += err * err;
            max_abs = max(max_abs, err);
            count++;
        }
    }

    for (int i = 0; i < N; ++i) {
        if (tr.theta_vertex_count[i] > 0) tr.theta_vertex[i] /= static_cast<double>(tr.theta_vertex_count[i]);
    }
    if (count > 0) {
        tr.theta_mean_abs = sum_abs / static_cast<double>(count);
        tr.theta_rms = sqrt(sum_sq / static_cast<double>(count));
        tr.theta_max_abs = max_abs;
    }
    return tr;
}

static map<int, int> degree_histogram(const vector<Face>& faces, int N) {
    auto edges = build_edges(faces);
    vector<int> degree(N, 0);
    for (const auto& kv : edges) {
        degree[kv.first.u]++;
        degree[kv.first.v]++;
    }
    map<int, int> hist;
    for (int d : degree) hist[d]++;
    return hist;
}

static bool finite_or_nan(double x) {
    return isfinite(x) || isnan(x);
}

static void write_outputs(const string& out_dir, const Config& cfg, int N, int flips_accepted,
                          const vector<Face>& faces, const vector<EdgeInfo>& edge_infos,
                          const RunResult& result, const EmbedResult& embed, const ThetaResult& theta,
                          const TriangulationStats& tri_stats) {
    filesystem::create_directories(out_dir);

    {
        ofstream f(filesystem::path(out_dir) / "stats.csv");
        f << "iter,E_K,max_abs_K,step,r_min,r_max,accepted\n";
        f << scientific << setprecision(17);
        for (const StatRow& row : result.stats) {
            f << row.iter << ',' << row.E_K << ',' << row.max_abs_K << ',' << row.step << ','
              << row.r_min << ',' << row.r_max << ',' << row.accepted << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "radii.csv");
        f << "vertex,u,r\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < N; ++i) f << i << ',' << result.u[i] << ',' << result.state.r[i] << '\n';
    }

    {
        ofstream f(filesystem::path(out_dir) / "faces.dat");
        for (const Face& face : faces) f << face.c[0].v << ' ' << face.c[1].v << ' ' << face.c[2].v << '\n';
    }

    {
        ofstream f(filesystem::path(out_dir) / "degree_hist.csv");
        f << "degree,count\n";
        for (const auto& kv : degree_histogram(faces, N)) f << kv.first << ',' << kv.second << '\n';
    }

    {
        ofstream f(filesystem::path(out_dir) / "edges.csv");
        f << "i,j,shift_x,shift_y,target_length\n";
        f << scientific << setprecision(17);
        for (const EdgeInfo& e : edge_infos) {
            f << e.i << ',' << e.j << ',' << e.shift_x << ',' << e.shift_y << ','
              << result.state.r[e.i] + result.state.r[e.j] << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "embed_stats.csv");
        f << "iter,E_embed,max_edge_abs_residual,rms_edge_residual,Lx,Ly,shear\n";
        f << scientific << setprecision(17);
        for (const EmbedStatsRow& row : embed.stats) {
            f << row.iter << ',' << row.E_embed << ',' << row.max_edge_abs_residual << ','
              << row.rms_edge_residual << ',' << row.Lx << ',' << row.Ly << ',' << row.shear << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "coords.csv");
        f << "vertex,sx,sy,x_real,y_real,u,r\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < N; ++i) {
            double sx = embed.ran ? embed.sx[i] : numeric_limits<double>::quiet_NaN();
            double sy = embed.ran ? embed.sy[i] : numeric_limits<double>::quiet_NaN();
            double x = embed.ran ? embed.Lx * sx + embed.shear * sy : numeric_limits<double>::quiet_NaN();
            double y = embed.ran ? embed.Ly * sy : numeric_limits<double>::quiet_NaN();
            f << i << ',' << sx << ',' << sy << ',' << x << ',' << y << ','
              << result.u[i] << ',' << result.state.r[i] << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "theta_vertex.csv");
        f << "vertex,Theta_i,incident_face_corners\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < N; ++i) {
            double tv = theta.theta_vertex.empty() ? numeric_limits<double>::quiet_NaN() : theta.theta_vertex[i];
            int cnt = theta.theta_vertex_count.empty() ? 0 : theta.theta_vertex_count[i];
            f << i << ',' << tv << ',' << cnt << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "theta_summary.txt");
        f << scientific << setprecision(17);
        f << "intrinsic_theta_residual: 0 by construction in the combinatorial tangent-triangle model\n";
        f << "coordinate_theta_mean_abs: " << theta.theta_mean_abs << "\n";
        f << "coordinate_theta_rms: " << theta.theta_rms << "\n";
        f << "coordinate_theta_max_abs: " << theta.theta_max_abs << "\n";
    }

    {
        bool strong_candidate = result.state.max_abs_K < cfg.tol &&
                                embed.ran && embed.max_edge_abs_residual < cfg.embed_tol &&
                                theta.theta_mean_abs < 1e-6;
        ofstream f(filesystem::path(out_dir) / "summary.txt");
        f << scientific << setprecision(17);
        f << "INC-Ricci v2 combinatorial plus coordinate proof-of-concept\n";
        f << "nx: " << cfg.nx << "\n";
        f << "ny: " << cfg.ny << "\n";
        f << "N: " << N << "\n";
        f << "F: " << tri_stats.F << "\n";
        f << "E: " << tri_stats.E << "\n";
        f << "Euler characteristic: " << tri_stats.chi << "\n";
        f << "flips requested: " << cfg.flips << "\n";
        f << "flips accepted: " << flips_accepted << "\n";
        f << "seed: " << cfg.seed << "\n";
        f << "method: " << cfg.method << "\n";
        f << "bounded mode: " << (cfg.bounded ? "yes" : "no") << "\n";
        f << "u_min bound: " << cfg.u_min << "\n";
        f << "u_max bound: " << cfg.u_max << "\n";
        f << "boundary rejections: " << result.boundary_rejections << "\n";
        f << "iterations: " << result.iterations << "\n";
        f << "initial E_K: " << result.initial_E_K << "\n";
        f << "initial max_abs_K: " << result.initial_max_abs_K << "\n";
        f << "final E_K: " << result.state.E_K << "\n";
        f << "final max_abs_K: " << result.state.max_abs_K << "\n";
        f << "final r_min: " << result.state.r_min << "\n";
        f << "final r_max: " << result.state.r_max << "\n";
        f << "final radius ratio: " << result.radius_ratio << "\n";
        f << "final u_min_actual: " << result.u_min_actual << "\n";
        f << "final u_max_actual: " << result.u_max_actual << "\n";
        f << "stopped reason: " << result.stopped_reason << "\n";
        f << "E_K monotonically decreased: " << (result.monotonic_E ? "yes" : "no") << "\n";
        f << "embed requested: " << (cfg.embed ? "yes" : "no") << "\n";
        f << "embed iterations: " << embed.iterations << "\n";
        f << "embed stopped reason: " << embed.stopped_reason << "\n";
        f << "final E_embed: " << embed.E_embed << "\n";
        f << "final max edge residual: " << embed.max_edge_abs_residual << "\n";
        f << "final rms edge residual: " << embed.rms_edge_residual << "\n";
        f << "final Lx: " << embed.Lx << "\n";
        f << "final Ly: " << embed.Ly << "\n";
        f << "final shear: " << embed.shear << "\n";
        f << "intrinsic_theta_residual: 0 by construction in the combinatorial tangent-triangle model\n";
        f << "coordinate theta_mean_abs: " << theta.theta_mean_abs << "\n";
        f << "coordinate theta_rms: " << theta.theta_rms << "\n";
        f << "coordinate theta_max_abs: " << theta.theta_max_abs << "\n";
        f << "strong_INC_candidate_thresholds: max_abs_K < tol, max_edge_abs_residual < embed_tol, theta_mean_abs < 1e-6\n";
        f << "strong_INC_candidate: " << (strong_candidate ? "yes" : "no") << "\n";
        f << "interpretation: intrinsic zero curvature alone is not a coordinate-level strong-INC success.\n";
        f << "limitation: This program tests a coordinate embedding by stress optimization but is not a full coordinate-level packing solver.\n";
        f << "limitation: It does not implement the Nature Materials SMC+FIRE protocol or prove physical convergence.\n";
    }
}

static void usage(const char* argv0) {
    cerr << "Usage: " << argv0
         << " [--nx N] [--ny N] [--flips N] [--seed N] [--max_iter N]"
         << " [--tol X] [--method ricci|grad] [--bounded 0|1] [--u_min X] [--u_max X]"
         << " [--embed 0|1] [--embed_iter N] [--embed_tol X] [--embed_lr X]"
         << " [--out DIR] [--test]\n";
}

static bool parse_bool_int(const string& s) {
    if (s == "1" || s == "true" || s == "yes") return true;
    if (s == "0" || s == "false" || s == "no") return false;
    throw runtime_error("expected boolean 0|1, got: " + s);
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        auto need_value = [&](const string& name) -> string {
            if (i + 1 >= argc) throw runtime_error("missing value for " + name);
            return argv[++i];
        };

        if (arg == "--nx") cfg.nx = stoi(need_value(arg));
        else if (arg == "--ny") cfg.ny = stoi(need_value(arg));
        else if (arg == "--flips") cfg.flips = stoi(need_value(arg));
        else if (arg == "--seed") cfg.seed = stoull(need_value(arg));
        else if (arg == "--max_iter") cfg.max_iter = stoi(need_value(arg));
        else if (arg == "--tol") cfg.tol = stod(need_value(arg));
        else if (arg == "--method") cfg.method = need_value(arg);
        else if (arg == "--bounded") cfg.bounded = parse_bool_int(need_value(arg));
        else if (arg == "--u_min") cfg.u_min = stod(need_value(arg));
        else if (arg == "--u_max") cfg.u_max = stod(need_value(arg));
        else if (arg == "--embed") cfg.embed = parse_bool_int(need_value(arg));
        else if (arg == "--embed_iter") cfg.embed_iter = stoi(need_value(arg));
        else if (arg == "--embed_tol") cfg.embed_tol = stod(need_value(arg));
        else if (arg == "--embed_lr") cfg.embed_lr = stod(need_value(arg));
        else if (arg == "--out") cfg.out = need_value(arg);
        else if (arg == "--test") cfg.test = true;
        else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            exit(0);
        } else {
            throw runtime_error("unknown argument: " + arg);
        }
    }

    if (cfg.nx <= 2 || cfg.ny <= 2) throw runtime_error("nx and ny must be greater than 2");
    if (cfg.flips < 0) throw runtime_error("flips must be non-negative");
    if (cfg.max_iter < 0) throw runtime_error("max_iter must be non-negative");
    if (cfg.tol <= 0.0) throw runtime_error("tol must be positive");
    if (cfg.method != "ricci" && cfg.method != "grad") throw runtime_error("method must be ricci or grad");
    if (cfg.u_min >= cfg.u_max) throw runtime_error("u_min must be less than u_max");
    if (cfg.embed_iter < 0) throw runtime_error("embed_iter must be non-negative");
    if (cfg.embed_tol <= 0.0) throw runtime_error("embed_tol must be positive");
    if (cfg.embed_lr <= 0.0) throw runtime_error("embed_lr must be positive");
    return cfg;
}

static RunResult run_method(const Config& cfg, vector<double> u, const vector<Face>& faces, bool verbose) {
    if (cfg.method == "ricci") {
        return run_ricci(std::move(u), faces, cfg.max_iter, cfg.tol, cfg.bounded, cfg.u_min, cfg.u_max, verbose);
    }
    if (cfg.method == "grad") return run_grad(std::move(u), faces, cfg.max_iter, cfg.tol, verbose);
    throw runtime_error("unknown method: " + cfg.method);
}

static bool run_tests() {
    bool ok = true;
    cout << "Running v2 tests\n";

    {
        int nx = 8, ny = 8;
        int N = nx * ny;
        auto faces = generate_periodic_triangular_torus(nx, ny);
        string err;
        TriangulationStats stats;
        bool valid = validate_triangulation(faces, N, &stats, &err);
        vector<double> u(N, 0.0);
        State s = compute_state(u, faces);
        string edge_err;
        vector<EdgeInfo> edges = build_edge_infos(faces, &edge_err);
        EmbedResult emb = run_embedding(edges, s.r, nx, ny, 80000, 1e-10, 0.02);
        ThetaResult th = compute_theta_residual(faces, s.r, emb);
        bool pass = valid && s.max_abs_K < 1e-12 && emb.max_edge_abs_residual < 1e-7 && th.theta_mean_abs < 1e-7;
        ok = ok && pass;
        cout << "Test 1 regular flat torus + embedding: " << (pass ? "PASS" : "FAIL")
             << " max_abs_K=" << scientific << setprecision(6) << s.max_abs_K
             << " max_edge_res=" << emb.max_edge_abs_residual
             << " theta_mean_abs=" << th.theta_mean_abs << '\n';
        if (!valid) cout << "  validation error: " << err << '\n';
        if (!edge_err.empty()) cout << "  edge error: " << edge_err << '\n';
    }

    RunResult bounded_result_for_embedding;
    vector<Face> bounded_faces_for_embedding;
    int bounded_nx = 8;
    int bounded_ny = 8;

    {
        int nx = 8, ny = 8;
        int N = nx * ny;
        mt19937_64 rng(1);
        auto faces = generate_periodic_triangular_torus(nx, ny);
        int accepted = 0;
        for (int i = 0; i < 200; ++i) {
            if (attempt_random_flip(faces, N, rng)) accepted++;
        }
        vector<double> u = random_log_radii(N, rng, 0.1);
        RunResult r = run_ricci(u, faces, 10000, 1e-8, false, -2.0, 2.0, false);
        bool pass = r.state.max_abs_K < 0.5 * r.initial_max_abs_K && r.state.E_K < r.initial_E_K && r.monotonic_E;
        ok = ok && pass;
        cout << "Test 2 randomized unbounded ricci: " << (pass ? "PASS" : "FAIL")
             << " flips_accepted=" << accepted
             << " initial_max_abs_K=" << r.initial_max_abs_K
             << " final_max_abs_K=" << r.state.max_abs_K << '\n';
    }

    {
        int nx = bounded_nx, ny = bounded_ny;
        int N = nx * ny;
        mt19937_64 rng(1);
        auto faces = generate_periodic_triangular_torus(nx, ny);
        int accepted = 0;
        for (int i = 0; i < 200; ++i) {
            if (attempt_random_flip(faces, N, rng)) accepted++;
        }
        vector<double> u = random_log_radii(N, rng, 0.1);
        RunResult r = run_ricci(u, faces, 10000, 1e-8, true, -2.0, 2.0, false);
        bool bounds_ok = r.u_min_actual >= -2.0 - 1e-12 && r.u_max_actual <= 2.0 + 1e-12;
        bool pass = bounds_ok && r.state.max_abs_K < 0.75 * r.initial_max_abs_K && r.state.E_K < r.initial_E_K;
        ok = ok && pass;
        bounded_result_for_embedding = r;
        bounded_faces_for_embedding = faces;
        cout << "Test 3 randomized bounded ricci: " << (pass ? "PASS" : "FAIL")
             << " flips_accepted=" << accepted
             << " final_max_abs_K=" << r.state.max_abs_K
             << " boundary_rejections=" << r.boundary_rejections
             << " u_range=[" << r.u_min_actual << "," << r.u_max_actual << "]"
             << " stopped=" << r.stopped_reason << '\n';
    }

    {
        string edge_err;
        vector<EdgeInfo> edges = build_edge_infos(bounded_faces_for_embedding, &edge_err);
        EmbedResult emb = run_embedding(edges, bounded_result_for_embedding.state.r, bounded_nx, bounded_ny,
                                        30000, 1e-8, 0.01);
        ThetaResult th = compute_theta_residual(bounded_faces_for_embedding, bounded_result_for_embedding.state.r, emb);
        bool pass = edge_err.empty() && finite_or_nan(emb.E_embed) && finite_or_nan(th.theta_mean_abs);
        ok = ok && pass;
        cout << "Test 4 bounded result embedding diagnostic: " << (pass ? "PASS" : "FAIL")
             << " E_embed=" << emb.E_embed
             << " max_edge_res=" << emb.max_edge_abs_residual
             << " theta_mean_abs=" << th.theta_mean_abs
             << " theta_max_abs=" << th.theta_max_abs << '\n';
    }

    cout << "Tests " << (ok ? "PASSED" : "FAILED") << '\n';
    return ok;
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);
        if (cfg.test) return run_tests() ? 0 : 1;

        const int N = cfg.nx * cfg.ny;
        mt19937_64 rng(cfg.seed);
        vector<Face> faces = generate_periodic_triangular_torus(cfg.nx, cfg.ny);

        string err;
        TriangulationStats tri_stats;
        if (!validate_triangulation(faces, N, &tri_stats, &err)) {
            throw runtime_error("initial triangulation invalid: " + err);
        }

        int flips_accepted = 0;
        for (int i = 0; i < cfg.flips; ++i) {
            if (attempt_random_flip(faces, N, rng)) flips_accepted++;
        }

        if (!validate_triangulation(faces, N, &tri_stats, &err)) {
            throw runtime_error("final triangulation invalid: " + err);
        }

        vector<double> u = random_log_radii(N, rng, 0.1);
        RunResult result = run_method(cfg, std::move(u), faces, true);

        string edge_err;
        vector<EdgeInfo> edge_infos = build_edge_infos(faces, &edge_err);
        if (!edge_err.empty()) throw runtime_error("edge shift validation failed: " + edge_err);

        EmbedResult embed;
        ThetaResult theta;
        if (cfg.embed) {
            cout << "embedding start edges " << edge_infos.size() << '\n';
            embed = run_embedding(edge_infos, result.state.r, cfg.nx, cfg.ny,
                                  cfg.embed_iter, cfg.embed_tol, cfg.embed_lr, true);
            theta = compute_theta_residual(faces, result.state.r, embed);
            cout << scientific << setprecision(17)
                 << "embedding_done E_embed " << embed.E_embed
                 << " max_edge_abs_residual " << embed.max_edge_abs_residual
                 << " theta_mean_abs " << theta.theta_mean_abs
                 << " theta_max_abs " << theta.theta_max_abs << '\n';
        }

        write_outputs(cfg.out, cfg, N, flips_accepted, faces, edge_infos, result, embed, theta, tri_stats);

        bool strong_candidate = result.state.max_abs_K < cfg.tol && cfg.embed &&
                                embed.max_edge_abs_residual < cfg.embed_tol &&
                                theta.theta_mean_abs < 1e-6;

        cout << scientific << setprecision(17);
        cout << "Done\n";
        cout << "out_dir " << cfg.out << '\n';
        cout << "flips_accepted " << flips_accepted << " / " << cfg.flips << '\n';
        cout << "final_E_K " << result.state.E_K << '\n';
        cout << "final_max_abs_K " << result.state.max_abs_K << '\n';
        cout << "radius_ratio " << result.radius_ratio << '\n';
        cout << "boundary_rejections " << result.boundary_rejections << '\n';
        cout << "final_max_edge_abs_residual " << embed.max_edge_abs_residual << '\n';
        cout << "theta_mean_abs " << theta.theta_mean_abs << '\n';
        cout << "theta_max_abs " << theta.theta_max_abs << '\n';
        cout << "strong_INC_candidate " << (strong_candidate ? "yes" : "no") << '\n';
        cout << "E_K_monotonic " << (result.monotonic_E ? "yes" : "no") << '\n';
        cout << "stopped_reason " << result.stopped_reason << '\n';
        cout << "embed_stopped_reason " << embed.stopped_reason << '\n';
        return 0;
    } catch (const exception& e) {
        cerr << "error: " << e.what() << '\n';
        usage(argv[0]);
        return 1;
    }
}
