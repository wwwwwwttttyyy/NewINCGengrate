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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

static constexpr double PI = 3.141592653589793238462643383279502884;

struct Face {
    int a, b, c;
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

struct RunResult {
    vector<double> u;
    State state;
    vector<StatRow> stats;
    int iterations = 0;
    string stopped_reason;
    bool monotonic_E = true;
    double initial_E_K = 0.0;
    double initial_max_abs_K = 0.0;
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

static vector<Face> generate_periodic_triangular_torus(int nx, int ny) {
    vector<Face> faces;
    faces.reserve(static_cast<size_t>(2 * nx * ny));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int A = vertex_id(i, j, nx, ny);
            int B = vertex_id(i + 1, j, nx, ny);
            int C = vertex_id(i, j + 1, nx, ny);
            int D = vertex_id(i + 1, j + 1, nx, ny);
            faces.push_back({A, B, C});
            faces.push_back({B, D, C});
        }
    }
    return faces;
}

static map<EdgeKey, vector<int>> build_edges(const vector<Face>& faces) {
    map<EdgeKey, vector<int>> edge_to_faces;
    for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi) {
        const Face& f = faces[fi];
        edge_to_faces[EdgeKey(f.a, f.b)].push_back(fi);
        edge_to_faces[EdgeKey(f.b, f.c)].push_back(fi);
        edge_to_faces[EdgeKey(f.c, f.a)].push_back(fi);
    }
    return edge_to_faces;
}

static bool validate_triangulation(const vector<Face>& faces, int V, TriangulationStats* out_stats = nullptr,
                                   string* error = nullptr) {
    auto fail = [&](const string& msg) {
        if (error) *error = msg;
        return false;
    };

    for (size_t i = 0; i < faces.size(); ++i) {
        const Face& f = faces[i];
        if (f.a < 0 || f.a >= V || f.b < 0 || f.b >= V || f.c < 0 || f.c >= V) {
            return fail("face has out-of-range vertex index");
        }
        if (f.a == f.b || f.b == f.c || f.c == f.a) {
            return fail("face has repeated vertices");
        }
    }

    auto edge_to_faces = build_edges(faces);
    for (const auto& kv : edge_to_faces) {
        if (kv.second.size() != 2) {
            return fail("edge does not have exactly two incident faces");
        }
    }

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

static int opposite_vertex(const Face& f, int a, int b) {
    if (f.a != a && f.a != b) return f.a;
    if (f.b != a && f.b != b) return f.b;
    if (f.c != a && f.c != b) return f.c;
    return -1;
}

static bool nondegenerate(const Face& f) {
    return f.a != f.b && f.b != f.c && f.c != f.a;
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

    int f1 = inc[0];
    int f2 = inc[1];
    int a = e.u;
    int b = e.v;
    int c = opposite_vertex(faces[f1], a, b);
    int d = opposite_vertex(faces[f2], a, b);
    if (c < 0 || d < 0 || c == d) return false;
    if (edge_to_faces.find(EdgeKey(c, d)) != edge_to_faces.end()) return false;

    Face nf1{c, d, a};
    Face nf2{d, c, b};
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
        if (scale <= 0.0 || !isfinite(scale)) {
            return numeric_limits<double>::quiet_NaN();
        }
        x /= scale;
        y /= scale;
        z /= scale;
        double denom = 2.0 * x * y;
        if (denom <= 0.0) return numeric_limits<double>::quiet_NaN();
        double c = (x * x + y * y - z * z) / denom;
        return acos(clamp_unit(c));
    };

    return {angle_at(ra, rb, rc), angle_at(rb, rc, ra), angle_at(rc, ra, rb)};
}

static array<double, 3> compute_angles_for_face(const Face& f, const vector<double>& r) {
    return compute_angles_from_radii(r[f.a], r[f.b], r[f.c]);
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
        s.K[f.a] -= ang[0];
        s.K[f.b] -= ang[1];
        s.K[f.c] -= ang[2];
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

static vector<LocalStencil> build_local_stencils(const vector<Face>& faces, int N) {
    vector<vector<int>> incident(N);
    for (int fi = 0; fi < static_cast<int>(faces.size()); ++fi) {
        const Face& f = faces[fi];
        incident[f.a].push_back(fi);
        incident[f.b].push_back(fi);
        incident[f.c].push_back(fi);
    }

    vector<LocalStencil> stencils(N);
    for (int i = 0; i < N; ++i) {
        vector<int> affected;
        for (int fi : incident[i]) {
            const Face& f = faces[fi];
            int verts[3] = {f.a, f.b, f.c};
            for (int v : verts) {
                if (find(affected.begin(), affected.end(), v) == affected.end()) {
                    affected.push_back(v);
                }
            }
        }
        stencils[i].affected_vertices = affected;
        for (int fi : incident[i]) {
            const Face& f = faces[fi];
            int verts[3] = {f.a, f.b, f.c};
            LocalFaceRef ref;
            ref.face = fi;
            for (int q = 0; q < 3; ++q) {
                auto it = find(affected.begin(), affected.end(), verts[q]);
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
            double rp[3] = {state.r[f.a], state.r[f.b], state.r[f.c]};
            double rm[3] = {state.r[f.a], state.r[f.b], state.r[f.c]};
            if (f.a == i) {
                rp[0] = ri_plus;
                rm[0] = ri_minus;
            } else if (f.b == i) {
                rp[1] = ri_plus;
                rm[1] = ri_minus;
            } else {
                rp[2] = ri_plus;
                rm[2] = ri_minus;
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

static RunResult run_ricci(vector<double> u, const vector<Face>& faces, int max_iter, double tol, bool verbose) {
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
                return result;
            }

            for (size_t i = 0; i < result.u.size(); ++i) {
                candidate_u[i] = result.u[i] - dt * state.K[i];
            }
            normalize_u(candidate_u);
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
            return result;
        }
    }

    result.state = state;
    result.stopped_reason = "reached max_iter";
    if (result.stats.empty() || result.stats.back().iter != result.iterations) {
        result.stats.push_back(make_stat_row(result.iterations, result.state, dt, 1));
    }
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
                return result;
            }

            for (size_t i = 0; i < result.u.size(); ++i) {
                candidate_u[i] = result.u[i] - trial_alpha * grad[i];
            }
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
            return result;
        }
    }

    result.state = state;
    result.stopped_reason = "reached max_iter";
    if (result.stats.empty() || result.stats.back().iter != result.iterations) {
        result.stats.push_back(make_stat_row(result.iterations, result.state, last_alpha, 1));
    }
    return result;
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

static void write_outputs(const string& out_dir, const Config& cfg, int N, int flips_accepted,
                          const vector<Face>& faces, const RunResult& result,
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
        for (int i = 0; i < N; ++i) {
            f << i << ',' << result.u[i] << ',' << result.state.r[i] << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "faces.dat");
        for (const Face& face : faces) {
            f << face.a << ' ' << face.b << ' ' << face.c << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "degree_hist.csv");
        f << "degree,count\n";
        for (const auto& kv : degree_histogram(faces, N)) {
            f << kv.first << ',' << kv.second << '\n';
        }
    }

    {
        ofstream f(filesystem::path(out_dir) / "summary.txt");
        f << scientific << setprecision(17);
        f << "INC-Ricci combinatorial proof-of-concept\n";
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
        f << "iterations: " << result.iterations << "\n";
        f << "initial E_K: " << result.initial_E_K << "\n";
        f << "initial max_abs_K: " << result.initial_max_abs_K << "\n";
        f << "final E_K: " << result.state.E_K << "\n";
        f << "final max_abs_K: " << result.state.max_abs_K << "\n";
        f << "final r_min: " << result.state.r_min << "\n";
        f << "final r_max: " << result.state.r_max << "\n";
        f << "stopped reason: " << result.stopped_reason << "\n";
        f << "E_K monotonically decreased: " << (result.monotonic_E ? "yes" : "no") << "\n";
        f << "theta_residual: 0 by construction in the combinatorial tangent-triangle model\n";
        f << "curvature_residual: " << result.state.max_abs_K << "\n";
        f << "interpretation: curvature_residual -> 0 means the tangent triangles close consistently into a flat torus circle-packing metric.\n";
        f << "strong INC endpoint: theta_residual = 0 and curvature_residual -> 0 in this combinatorial model.\n";
        f << "limitation: This program tests the combinatorial radius/curvature part of the INC-Ricci idea.\n";
        f << "limitation: It does not reconstruct explicit Euclidean coordinates of circle centers.\n";
        f << "limitation: It does not prove that the Nature Materials SMC+FIRE protocol converges.\n";
        f << "limitation: It tests whether a randomized triangulated torus can be driven toward a zero-curvature tangent-circle metric.\n";
    }
}

static void usage(const char* argv0) {
    cerr << "Usage: " << argv0
         << " [--nx N] [--ny N] [--flips N] [--seed N] [--max_iter N]"
         << " [--tol X] [--method ricci|grad] [--out DIR] [--test]\n";
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        auto need_value = [&](const string& name) -> string {
            if (i + 1 >= argc) {
                throw runtime_error("missing value for " + name);
            }
            return argv[++i];
        };

        if (arg == "--nx") {
            cfg.nx = stoi(need_value(arg));
        } else if (arg == "--ny") {
            cfg.ny = stoi(need_value(arg));
        } else if (arg == "--flips") {
            cfg.flips = stoi(need_value(arg));
        } else if (arg == "--seed") {
            cfg.seed = stoull(need_value(arg));
        } else if (arg == "--max_iter") {
            cfg.max_iter = stoi(need_value(arg));
        } else if (arg == "--tol") {
            cfg.tol = stod(need_value(arg));
        } else if (arg == "--method") {
            cfg.method = need_value(arg);
        } else if (arg == "--out") {
            cfg.out = need_value(arg);
        } else if (arg == "--test") {
            cfg.test = true;
        } else if (arg == "--help" || arg == "-h") {
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
    return cfg;
}

static RunResult run_method(const string& method, vector<double> u, const vector<Face>& faces,
                            int max_iter, double tol, bool verbose) {
    if (method == "ricci") return run_ricci(std::move(u), faces, max_iter, tol, verbose);
    if (method == "grad") return run_grad(std::move(u), faces, max_iter, tol, verbose);
    throw runtime_error("unknown method: " + method);
}

static bool run_tests() {
    bool ok = true;
    cout << "Running tests\n";

    {
        int nx = 8, ny = 8;
        int N = nx * ny;
        auto faces = generate_periodic_triangular_torus(nx, ny);
        string err;
        TriangulationStats stats;
        bool valid = validate_triangulation(faces, N, &stats, &err);
        vector<double> u(N, 0.0);
        State s = compute_state(u, faces);
        bool pass = valid && s.max_abs_K < 1e-12;
        ok = ok && pass;
        cout << "Test 1 regular flat torus: " << (pass ? "PASS" : "FAIL")
             << " max_abs_K=" << scientific << setprecision(6) << s.max_abs_K << '\n';
        if (!valid) cout << "  validation error: " << err << '\n';
    }

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
        RunResult r = run_ricci(u, faces, 10000, 1e-8, false);
        bool pass = r.state.max_abs_K < 0.5 * r.initial_max_abs_K && r.state.E_K < r.initial_E_K && r.monotonic_E;
        ok = ok && pass;
        cout << "Test 2 randomized ricci: " << (pass ? "PASS" : "FAIL")
             << " flips_accepted=" << accepted
             << " initial_max_abs_K=" << r.initial_max_abs_K
             << " final_max_abs_K=" << r.state.max_abs_K
             << " initial_E_K=" << r.initial_E_K
             << " final_E_K=" << r.state.E_K << '\n';
    }

    {
        int nx = 8, ny = 8;
        int N = nx * ny;
        mt19937_64 rng(1);
        auto faces = generate_periodic_triangular_torus(nx, ny);
        for (int i = 0; i < 200; ++i) {
            attempt_random_flip(faces, N, rng);
        }
        vector<double> u = random_log_radii(N, rng, 0.1);
        RunResult r = run_grad(u, faces, 500, 1e-8, false);
        bool pass = r.monotonic_E && r.state.E_K < r.initial_E_K;
        ok = ok && pass;
        cout << "Test 3 randomized grad monotonic: " << (pass ? "PASS" : "FAIL")
             << " monotonic_E=" << (r.monotonic_E ? "yes" : "no")
             << " initial_E_K=" << r.initial_E_K
             << " final_E_K=" << r.state.E_K
             << " final_max_abs_K=" << r.state.max_abs_K << '\n';
    }

    cout << "Tests " << (ok ? "PASSED" : "FAILED") << '\n';
    return ok;
}

int main(int argc, char** argv) {
    try {
        Config cfg = parse_args(argc, argv);
        if (cfg.test) {
            return run_tests() ? 0 : 1;
        }

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
        RunResult result = run_method(cfg.method, std::move(u), faces, cfg.max_iter, cfg.tol, true);

        write_outputs(cfg.out, cfg, N, flips_accepted, faces, result, tri_stats);

        cout << scientific << setprecision(17);
        cout << "Done\n";
        cout << "out_dir " << cfg.out << '\n';
        cout << "flips_accepted " << flips_accepted << " / " << cfg.flips << '\n';
        cout << "final_E_K " << result.state.E_K << '\n';
        cout << "final_max_abs_K " << result.state.max_abs_K << '\n';
        cout << "E_K_monotonic " << (result.monotonic_E ? "yes" : "no") << '\n';
        cout << "stopped_reason " << result.stopped_reason << '\n';
        return 0;
    } catch (const exception& e) {
        cerr << "error: " << e.what() << '\n';
        usage(argv[0]);
        return 1;
    }
}
