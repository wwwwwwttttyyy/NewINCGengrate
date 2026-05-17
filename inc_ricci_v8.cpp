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
    double rms_K = 0.0;
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
    double final_E_reg = numeric_limits<double>::quiet_NaN();
};

struct DistributionDiagnostics {
    bool ran = false;
    vector<double> u0;
    vector<double> r0;
    double D_u = 0.0;
    double D_sort = 0.0;
    double D_quantile = 0.0;
    double D_delta = 0.0;
    double D_ratio = 0.0;
    double q05_u = numeric_limits<double>::quiet_NaN();
    double q50_u = numeric_limits<double>::quiet_NaN();
    double q95_u = numeric_limits<double>::quiet_NaN();
    double q05_u0 = numeric_limits<double>::quiet_NaN();
    double q50_u0 = numeric_limits<double>::quiet_NaN();
    double q95_u0 = numeric_limits<double>::quiet_NaN();
    double current_delta = numeric_limits<double>::quiet_NaN();
    double target_delta = numeric_limits<double>::quiet_NaN();
    double current_radius_ratio = numeric_limits<double>::quiet_NaN();
    double target_radius_ratio = numeric_limits<double>::quiet_NaN();
    bool fixed_intrinsic_good = false;
    bool weak_distribution_preserved = false;
    bool distribution_preserved = false;
};

struct LinearResponseResult {
    bool ran = false;
    string status = "not requested";
    vector<double> K0;
    vector<double> K_lin;
    vector<double> du;
    vector<double> vertex_mismatch;
    vector<double> edge_mismatch;
    double rms_K0 = numeric_limits<double>::quiet_NaN();
    double max_K0 = numeric_limits<double>::quiet_NaN();
    double rms_K_lin = numeric_limits<double>::quiet_NaN();
    double max_K_lin = numeric_limits<double>::quiet_NaN();
    double du_rms = numeric_limits<double>::quiet_NaN();
    double du_max = numeric_limits<double>::quiet_NaN();
    double du_l2 = numeric_limits<double>::quiet_NaN();
    double response_ratio = numeric_limits<double>::quiet_NaN();
    double S_pred = numeric_limits<double>::quiet_NaN();
    double degree_radius_corr = numeric_limits<double>::quiet_NaN();
    bool singular_response_flag = false;
    bool linear_predicted_good = false;
    int cg_iters = 0;
    double cg_residual = numeric_limits<double>::quiet_NaN();
};

struct DegreeDiagnostics {
    bool ran = false;
    vector<int> target_degree;
    double target_degree_radius_corr = numeric_limits<double>::quiet_NaN();
    double S_degree = numeric_limits<double>::quiet_NaN();
    double E_deg = numeric_limits<double>::quiet_NaN();
    int max_degree_error = 0;
    double degree_match_fraction = numeric_limits<double>::quiet_NaN();
    double degree_radius_corr = numeric_limits<double>::quiet_NaN();
    map<int, int> target_hist;
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
    int face_degenerate_count = 0;
    int face_negative_signed_count = 0;
    int orientation_neighbor_violation_count = 0;

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
    bool fixed_intrinsic_good = false;
    bool weak_distribution_preserved = false;
    bool distribution_preserved = false;
    bool physical_INC_like_candidate_v6 = false;

    double degree6_fraction = numeric_limits<double>::quiet_NaN();
    double edge_length_cv = numeric_limits<double>::quiet_NaN();
    bool noncrystalline_topology_heuristic = false;
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
    bool require_consistent_orientation = false;
    bool check_orientation_neighbors = true;
    bool check_crossings = true;
    bool check_crossings_explicit = false;
    double intersect_tol = 1e-10;
    int hist_bins = 50;
    double target_max_ratio = 10.0;
    double target_delta = 0.25;
    string radius_mode = "free";
    string target_dist = "lognormal";
    double target_radius_ratio = 5.0;
    unsigned long long radius_seed = 123;
    double radius_noise_std = 0.25;
    double fixed_tol_K = 1e-3;
    double fixed_tol_rms_K = 1e-4;
    double lambda_u = 0.0;
    double lambda_sort = 100.0;
    double lambda_quantile = 10.0;
    double lambda_delta = 100.0;
    double lambda_ratio = 100.0;
    double dist_tol_sort = 1e-4;
    double dist_tol_delta = 0.02;
    double dist_tol_ratio = 0.05;
    double develop_spread_tol = 1e-6;
    double period_fit_tol = 1e-6;
    double theta_tol = 1e-8;
    double contact_tol = 1e-6;

    bool surgery = false;
    int surgery_steps = 1000;
    int surgery_trials_per_step = 20;
    string surgery_metric = "v6";
    double surgery_temperature = 0.0;
    int surgery_ricci_iter = 20000;
    int surgery_report_every = 10;
    int surgery_batch_flips = 1;
    string proposal_mode = "random";

    double K_scale = 1e-8;
    double rmsK_scale = 1e-9;
    double dev_scale = 1e-6;
    double period_scale = 1e-6;
    double overlap_scale = 1e-8;
    double wK = 1.0;
    double wRmsK = 1.0;
    double wDev = 1.0;
    double wPer = 1.0;
    double wOv = 5.0;
    double wOvN = 5.0;
    double wCross = 10.0;
    double wDeg = 10.0;
    double wOrient = 3.0;
    double wRad = 5.0;
    double wDelta = 50.0;
    double wSort = 200.0;
    double wRatio = 50.0;

    bool multistart_surgery = false;
    int multistart_count = 10;
    string multistart_out = "multistart.csv";

    bool pareto_scan = false;
    string pareto_out = "pareto.csv";
    vector<double> pareto_delta_list;
    vector<double> pareto_ratio_list;
    vector<string> pareto_modes;
    vector<unsigned long long> pareto_seeds;

    double lin_h = 1e-6;
    double lin_lambda = 1e-2;
    double lin_cg_tol = 1e-10;
    int lin_cg_max_iter = 5000;
    double du_max_allowed = 0.1;
    double du_rms_allowed = 0.03;
    double K_rms_scale = 1e-3;
    double K_max_scale = 1e-2;
    double w_res = 1.0;
    double w_max = 1.0;
    double w_du = 2.0;
    double w_dumax = 2.0;
    double w_sing = 10.0;
    double w_mis = 1.0;
    double corr_target = 0.6;
    double mismatch_a = 1.0;
    double mismatch_b = 1.0;
    double mismatch_c = 0.5;

    string degree_target_mode = "none";
    double frac_deg4 = 0.08;
    double frac_deg5 = 0.17;
    double frac_deg6 = 0.50;
    double frac_deg7 = 0.17;
    double frac_deg8 = 0.08;
    double w_deg = 10.0;
    double w_degmax = 1.0;
    double w_corr = 10.0;
    bool degree_precondition = false;
    int degree_pre_steps = 1000;
    int degree_pre_trials = 50;
    double degree_pre_temperature = 0.0;
    string degree_proposal_mode = "greedy";
    double alpha_degree = 1.0;
    bool linear_surgery_use_degree = true;
    double beta_degree = 5.0;

    bool degree_fraction_scan = false;
    string degree_fraction_out = "degree_fraction_scan.csv";
    bool compare_topology_strategies = false;
    string strategy_out = "strategy_compare.csv";

    int linear_surgery_steps = 500;
    int linear_surgery_trials_per_step = 50;
    double linear_surgery_temperature = 0.0;

    bool refine_after_linear = false;
    int refine_max_iter = 50000;
    double refine_tol = 1e-8;
    double refine_du_bound = 0.15;
    double predictor_rms_tol = 1e-3;

    bool screen_topologies = false;
    int screen_count = 1000;
    int screen_flips = 100;
    string screen_out = "topology_screen.csv";
    bool screen_same_radii = true;

    bool linear_pareto_scan = false;
};

struct RunOutput {
    Config cfg;
    Mesh mesh;
    RicciResult ricci;
    DevelopResult develop;
    PackingDiagnostics packing;
    DistributionDiagnostics dist;
    LinearResponseResult linear;
    DegreeDiagnostics degree;
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
    int face_degenerate_count = 0;
    int orientation_neighbor_violation_count = 0;
    double polydispersity_delta = numeric_limits<double>::quiet_NaN();
    bool physically_reasonable_radius = false;
    bool valid_disk_packing_candidate = false;
    bool strong_INC_candidate = false;
    bool physical_INC_like_candidate = false;
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
    int face_degenerate_count = 0;
    int orientation_neighbor_violation_count = 0;
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

static bool attempt_flip_edge_index(Mesh& mesh, int edge_idx, const unordered_map<uint64_t, int>& edge_index) {
    if (edge_idx < 0 || edge_idx >= static_cast<int>(mesh.edges.size())) return false;
    const EdgeAdj edge = mesh.edges[edge_idx];
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

static bool attempt_random_flip(Mesh& mesh, mt19937_64& rng) {
    unordered_map<uint64_t, int> edge_index;
    string error;
    if (!build_adjacency(mesh, &error, &edge_index)) return false;
    if (mesh.edges.empty()) return false;

    uniform_int_distribution<int> pick(0, static_cast<int>(mesh.edges.size()) - 1);
    return attempt_flip_edge_index(mesh, pick(rng), edge_index);
}

static vector<int> vertex_degrees(const Mesh& mesh);

static bool attempt_biased_flip(Mesh& mesh, mt19937_64& rng, const Config& cfg, const RunOutput& current) {
    if (cfg.proposal_mode == "random") return attempt_random_flip(mesh, rng);
    unordered_map<uint64_t, int> edge_index;
    string error;
    if (!build_adjacency(mesh, &error, &edge_index)) return false;
    if (mesh.edges.empty()) return false;

    vector<int> candidates;
    auto add_edges_near_vertex = [&](int v) {
        for (int ei = 0; ei < static_cast<int>(mesh.edges.size()); ++ei) {
            if (mesh.edges[ei].a == v || mesh.edges[ei].b == v) candidates.push_back(ei);
        }
    };

    if (cfg.proposal_mode == "overlap" && !current.packing.overlap_samples.empty()) {
        const OverlapSample& s = current.packing.overlap_samples.front();
        add_edges_near_vertex(s.i);
        add_edges_near_vertex(s.j);
    }

    if (candidates.empty()) {
        const vector<double>& K = current.ricci.state.K;
        if (!K.empty()) {
            double total = 0.0;
            for (double k : K) total += fabs(k) + 1e-12;
            uniform_real_distribution<double> uni(0.0, total);
            double t = uni(rng);
            int picked_v = 0;
            for (int i = 0; i < static_cast<int>(K.size()); ++i) {
                t -= fabs(K[i]) + 1e-12;
                if (t <= 0.0) {
                    picked_v = i;
                    break;
                }
            }
            add_edges_near_vertex(picked_v);
        }
    }

    if (candidates.empty()) return attempt_random_flip(mesh, rng);
    shuffle(candidates.begin(), candidates.end(), rng);
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    for (int ei : candidates) {
        if (attempt_flip_edge_index(mesh, ei, edge_index)) return true;
    }
    return attempt_random_flip(mesh, rng);
}

static bool attempt_linear_guided_flip(Mesh& mesh, mt19937_64& rng, const Config& cfg,
                                       const RunOutput& current, int& flipped_i, int& flipped_j) {
    flipped_i = -1;
    flipped_j = -1;
    unordered_map<uint64_t, int> edge_index;
    string error;
    if (!build_adjacency(mesh, &error, &edge_index)) return false;
    if (mesh.edges.empty()) return false;

    vector<double> weights(mesh.edges.size(), 1.0);
    if (cfg.proposal_mode == "mismatch" && current.linear.edge_mismatch.size() == mesh.edges.size()) {
        for (int i = 0; i < static_cast<int>(mesh.edges.size()); ++i) {
            weights[i] = max(0.0, current.linear.edge_mismatch[i]) + 1e-12;
        }
    } else if (cfg.proposal_mode == "highK" && !current.linear.K_lin.empty()) {
        for (int i = 0; i < static_cast<int>(mesh.edges.size()); ++i) {
            const EdgeAdj& e = mesh.edges[i];
            weights[i] = fabs(current.linear.K_lin[e.a]) + fabs(current.linear.K_lin[e.b]) + 1e-12;
        }
    }
    if (cfg.linear_surgery_use_degree && current.degree.ran &&
        static_cast<int>(current.degree.target_degree.size()) == mesh.N) {
        vector<int> degree = vertex_degrees(mesh);
        for (int i = 0; i < static_cast<int>(mesh.edges.size()); ++i) {
            const EdgeAdj& e = mesh.edges[i];
            int c = opposite_index(mesh.faces[e.f0], e.a, e.b);
            int d = opposite_index(mesh.faces[e.f1], e.a, e.b);
            if (c < 0 || d < 0) continue;
            int vc = mesh.faces[e.f0].c[c].v;
            int vd = mesh.faces[e.f1].c[d].v;
            auto sqerr = [&](int v, int delta) {
                double before = static_cast<double>(degree[v] - current.degree.target_degree[v]);
                double after = static_cast<double>(degree[v] + delta - current.degree.target_degree[v]);
                return before * before - after * after;
            };
            double potential = sqerr(e.a, -1) + sqerr(e.b, -1) + sqerr(vc, +1) + sqerr(vd, +1);
            weights[i] += cfg.beta_degree * max(0.0, potential);
        }
    }

    double total = accumulate(weights.begin(), weights.end(), 0.0);
    if (!(total > 0.0) || !isfinite(total)) return attempt_random_flip(mesh, rng);
    uniform_real_distribution<double> uni(0.0, total);
    for (int attempt = 0; attempt < 40; ++attempt) {
        double t = uni(rng);
        int chosen = static_cast<int>(mesh.edges.size()) - 1;
        for (int i = 0; i < static_cast<int>(weights.size()); ++i) {
            t -= weights[i];
            if (t <= 0.0) {
                chosen = i;
                break;
            }
        }
        EdgeAdj e = mesh.edges[chosen];
        if (attempt_flip_edge_index(mesh, chosen, edge_index)) {
            flipped_i = e.a;
            flipped_j = e.b;
            return true;
        }
    }
    return false;
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
    s.rms_K = N > 0 ? sqrt((2.0 * s.E_K) / static_cast<double>(N)) : 0.0;
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

static DistributionDiagnostics compute_distribution_diagnostics(const vector<double>& u,
                                                                const vector<double>& u0,
                                                                const Config& cfg,
                                                                double max_abs_K,
                                                                double rms_K);
static double regularized_energy(const State& s, const vector<double>& u0, const Config& cfg);

static void finalize_ricci(RicciResult& rr) {
    if (!rr.state.u.empty()) {
        auto mm = minmax_element(rr.state.u.begin(), rr.state.u.end());
        rr.u_min_actual = *mm.first;
        rr.u_max_actual = *mm.second;
    }
    rr.radius_ratio = rr.state.r_min > 0.0 ? rr.state.r_max / rr.state.r_min
                                            : numeric_limits<double>::infinity();
}

static RicciResult run_ricci(const Mesh& mesh, vector<double> u, const Config& cfg,
                             const vector<double>& u0, bool keep_stats, bool verbose) {
    normalize_u(u);
    RicciResult rr;
    State current, candidate;
    compute_state_inplace(mesh, u, current);
    rr.initial_E_K = current.E_K;
    rr.initial_max_abs_K = current.max_abs_K;
    double current_E_reg = regularized_energy(current, u0, cfg);
    rr.final_E_reg = current_E_reg;
    if (keep_stats) rr.stats.push_back(make_stat(0, current, 0.01, 1));

    if (cfg.radius_mode == "fixed") {
        rr.state = std::move(current);
        rr.iterations = 0;
        rr.stopped_reason = "fixed radii: curvature evaluated without Ricci update";
        finalize_ricci(rr);
        rr.final_E_reg = rr.state.E_K;
        return rr;
    }

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

    auto weak_preserved_now = [&]() {
        if (cfg.radius_mode != "weak") return true;
        DistributionDiagnostics dd = compute_distribution_diagnostics(current.u, u0, cfg, current.max_abs_K, current.rms_K);
        return dd.weak_distribution_preserved;
    };

    if (current.max_abs_K < cfg.tol && weak_preserved_now()) {
        rr.state = std::move(current);
        rr.stopped_reason = "converged: max_abs_K < tol at iter 0";
        finalize_ricci(rr);
        rr.final_E_reg = current_E_reg;
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
                rr.final_E_reg = current_E_reg;
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
                rr.final_E_reg = current_E_reg;
                return rr;
            }

            double candidate_E_reg = regularized_energy(candidate, u0, cfg);
            if (candidate_E_reg < current_E_reg) {
                accepted = true;
                accepted_step = dt;
                rr.monotonic_E = rr.monotonic_E && (candidate.E_K <= current.E_K + 1e-15);
                u.swap(cand_u);
                current = std::move(candidate);
                current_E_reg = candidate_E_reg;
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

        if (current.max_abs_K < cfg.tol && weak_preserved_now()) {
            rr.state = std::move(current);
            rr.stopped_reason = cfg.radius_mode == "weak"
                                  ? "converged: max_abs_K < tol and distribution preserved"
                                  : "converged: max_abs_K < tol";
            finalize_ricci(rr);
            rr.final_E_reg = current_E_reg;
            return rr;
        }

        if (cfg.stagnation_window > 0 && iter - window_start >= cfg.stagnation_window) {
            double rel = (window_E - current.E_K) / max(window_E, 1e-300);
            if (rel < cfg.stagnation_rel) {
                rr.state = std::move(current);
                rr.stopped_reason = "stalled_by_stagnation";
                finalize_ricci(rr);
                rr.final_E_reg = current_E_reg;
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
    rr.final_E_reg = current_E_reg;
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

static vector<double> radii_from_u(const vector<double>& u) {
    vector<double> r(u.size());
    for (size_t i = 0; i < u.size(); ++i) r[i] = safe_exp(u[i]);
    return r;
}

static double radius_ratio_from_u(const vector<double>& u) {
    vector<double> r = radii_from_u(u);
    if (r.empty()) return numeric_limits<double>::quiet_NaN();
    auto mm = minmax_element(r.begin(), r.end());
    return *mm.first > 0.0 ? *mm.second / *mm.first : numeric_limits<double>::infinity();
}

static vector<double> generate_target_log_radii(int N, const Config& cfg) {
    vector<double> best;
    mt19937_64 rng(cfg.radius_seed);
    double max_ratio = cfg.target_radius_ratio > 0.0 ? cfg.target_radius_ratio : numeric_limits<double>::infinity();
    int attempts = isfinite(max_ratio) ? 1000 : 1;
    for (int attempt = 0; attempt < attempts; ++attempt) {
        vector<double> u(N, 0.0);
        if (cfg.target_dist == "lognormal") {
            double sigma = cfg.target_delta > 0.0 ? sqrt(log(1.0 + cfg.target_delta * cfg.target_delta)) : 0.0;
            normal_distribution<double> normal(0.0, sigma);
            for (double& x : u) x = normal(rng);
        } else if (cfg.target_dist == "uniform_u") {
            double B = 0.5 * log(max(cfg.target_radius_ratio, 1.0));
            uniform_real_distribution<double> uni(-B, B);
            for (double& x : u) x = uni(rng);
        } else if (cfg.target_dist == "powerlaw_like") {
            double R = max(cfg.target_radius_ratio, 1.0 + 1e-12);
            double rmin = 1.0;
            double rmax = R;
            uniform_real_distribution<double> uni(0.0, 1.0);
            for (double& x : u) {
                double z = uni(rng);
                double inv = (1.0 / rmin) - z * ((1.0 / rmin) - (1.0 / rmax));
                double r = 1.0 / max(inv, 1e-300);
                x = log(r);
            }
        } else {
            throw runtime_error("target_dist must be lognormal, uniform_u, or powerlaw_like");
        }
        normalize_u(u);
        double ratio = radius_ratio_from_u(u);
        if (!isfinite(max_ratio) || ratio <= max_ratio * (1.0 + 1e-12)) return u;
        if (best.empty() || ratio < radius_ratio_from_u(best)) best = u;
    }
    throw runtime_error("failed to sample target radii within target_radius_ratio after 1000 attempts");
}

static DistributionDiagnostics compute_distribution_diagnostics(const vector<double>& u,
                                                                const vector<double>& u0,
                                                                const Config& cfg,
                                                                double max_abs_K,
                                                                double rms_K) {
    DistributionDiagnostics dd;
    dd.ran = true;
    dd.u0 = u0;
    dd.r0 = radii_from_u(u0);
    dd.target_delta = cfg.target_delta;
    dd.target_radius_ratio = cfg.target_radius_ratio;
    if (u.empty()) return dd;

    vector<double> r = radii_from_u(u);
    RadiusStats cur = compute_radius_stats(r, u, cfg.hist_bins);
    RadiusStats tar = compute_radius_stats(dd.r0, u0, cfg.hist_bins);
    dd.current_delta = cur.polydispersity_delta;
    dd.current_radius_ratio = radius_ratio_from_u(u);
    dd.q05_u = cur.q05_u;
    dd.q50_u = cur.q50_u;
    dd.q95_u = cur.q95_u;
    dd.q05_u0 = tar.q05_u;
    dd.q50_u0 = tar.q50_u;
    dd.q95_u0 = tar.q95_u;

    if (!u0.empty() && u0.size() == u.size()) {
        double sum_du = 0.0;
        for (size_t i = 0; i < u.size(); ++i) {
            double d = u[i] - u0[i];
            sum_du += d * d;
        }
        dd.D_u = sum_du / static_cast<double>(u.size());

        vector<double> su = u;
        vector<double> su0 = u0;
        sort(su.begin(), su.end());
        sort(su0.begin(), su0.end());
        double sum_sort = 0.0;
        for (size_t i = 0; i < su.size(); ++i) {
            double d = su[i] - su0[i];
            sum_sort += d * d;
        }
        dd.D_sort = sum_sort / static_cast<double>(su.size());
        dd.D_quantile = pow(dd.q05_u - dd.q05_u0, 2.0) +
                        pow(dd.q50_u - dd.q50_u0, 2.0) +
                        pow(dd.q95_u - dd.q95_u0, 2.0);
    }
    dd.D_delta = pow(dd.current_delta - cfg.target_delta, 2.0);
    dd.D_ratio = pow(max(0.0, log(max(dd.current_radius_ratio / cfg.target_radius_ratio, 1e-300))), 2.0);

    dd.fixed_intrinsic_good = cfg.radius_mode == "fixed" &&
                              (max_abs_K < cfg.fixed_tol_K || rms_K < cfg.fixed_tol_rms_K);
    dd.weak_distribution_preserved = cfg.radius_mode == "weak" &&
                                     dd.D_sort < cfg.dist_tol_sort &&
                                     fabs(dd.current_delta - cfg.target_delta) < cfg.dist_tol_delta &&
                                     dd.current_radius_ratio <= cfg.target_radius_ratio * (1.0 + cfg.dist_tol_ratio);
    dd.distribution_preserved = cfg.radius_mode == "fixed" || dd.weak_distribution_preserved;
    return dd;
}

static double regularized_energy(const State& s, const vector<double>& u0, const Config& cfg) {
    if (cfg.radius_mode != "weak") return s.E_K;
    DistributionDiagnostics dd = compute_distribution_diagnostics(s.u, u0, cfg, s.max_abs_K, s.rms_K);
    return s.E_K +
           cfg.lambda_u * dd.D_u +
           cfg.lambda_sort * dd.D_sort +
           cfg.lambda_quantile * dd.D_quantile +
           cfg.lambda_delta * dd.D_delta +
           cfg.lambda_ratio * dd.D_ratio;
}

static vector<int> vertex_degrees(const Mesh& mesh) {
    vector<int> d(mesh.N, 0);
    for (const EdgeAdj& e : mesh.edges) {
        d[e.a]++;
        d[e.b]++;
    }
    return d;
}

static double pearson_corr_degree_u(const vector<int>& degree, const vector<double>& u) {
    int N = static_cast<int>(u.size());
    if (N == 0) return numeric_limits<double>::quiet_NaN();
    double mean_d = 0.0, mean_u = 0.0;
    for (int i = 0; i < N; ++i) {
        mean_d += degree[i];
        mean_u += u[i];
    }
    mean_d /= N;
    mean_u /= N;
    double cov = 0.0, vard = 0.0, varu = 0.0;
    for (int i = 0; i < N; ++i) {
        double dd = degree[i] - mean_d;
        double du = u[i] - mean_u;
        cov += dd * du;
        vard += dd * dd;
        varu += du * du;
    }
    if (vard < 1e-30 && varu < 1e-30) return 1.0;
    if (vard < 1e-30 || varu < 1e-30) return 0.0;
    return cov / sqrt(vard * varu);
}

static void adjust_target_degree_sum(vector<int>& target, const vector<double>& u0, int desired_sum) {
    int sum = accumulate(target.begin(), target.end(), 0);
    while (sum < desired_sum) {
        int best = -1;
        for (int i = 0; i < static_cast<int>(target.size()); ++i) {
            if (target[i] >= 8) continue;
            if (best < 0 || u0[i] > u0[best]) best = i;
        }
        if (best < 0) break;
        target[best]++;
        sum++;
    }
    while (sum > desired_sum) {
        int best = -1;
        for (int i = 0; i < static_cast<int>(target.size()); ++i) {
            if (target[i] <= 4) continue;
            if (best < 0 || u0[i] < u0[best]) best = i;
        }
        if (best < 0) break;
        target[best]--;
        sum--;
    }
}

static vector<int> assign_target_degrees(const vector<double>& u0, const Config& cfg) {
    const int N = static_cast<int>(u0.size());
    vector<int> target(N, 6);
    if (cfg.degree_target_mode == "none" || N == 0) return target;

    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        if (u0[a] != u0[b]) return u0[a] < u0[b];
        return a < b;
    });

    if (cfg.degree_target_mode == "quantile") {
        array<double, 5> frac = {cfg.frac_deg4, cfg.frac_deg5, cfg.frac_deg6, cfg.frac_deg7, cfg.frac_deg8};
        array<int, 5> deg = {4, 5, 6, 7, 8};
        array<int, 5> count = {};
        array<double, 5> rem = {};
        int used = 0;
        for (int i = 0; i < 5; ++i) {
            double exact = frac[i] * N;
            count[i] = static_cast<int>(floor(exact));
            rem[i] = exact - count[i];
            used += count[i];
        }
        while (used < N) {
            int best = 0;
            for (int i = 1; i < 5; ++i) if (rem[i] > rem[best]) best = i;
            count[best]++;
            rem[best] = -1.0;
            used++;
        }
        while (used > N) {
            int best = 0;
            for (int i = 1; i < 5; ++i) if (count[i] > count[best]) best = i;
            count[best]--;
            used--;
        }
        int pos = 0;
        for (int k = 0; k < 5; ++k) {
            for (int c = 0; c < count[k] && pos < N; ++c, ++pos) target[order[pos]] = deg[k];
        }
    } else if (cfg.degree_target_mode == "linear") {
        double mean = accumulate(u0.begin(), u0.end(), 0.0) / max(1, N);
        double var = 0.0;
        for (double u : u0) var += (u - mean) * (u - mean);
        double sd = sqrt(max(var / max(1, N), 1e-30));
        for (int i = 0; i < N; ++i) {
            int d = static_cast<int>(llround(6.0 + (u0[i] - mean) / sd));
            target[i] = max(4, min(8, d));
        }
    } else {
        throw runtime_error("degree_target_mode must be quantile, linear, or none");
    }
    adjust_target_degree_sum(target, u0, 6 * N);
    return target;
}

static DegreeDiagnostics compute_degree_diagnostics(const Mesh& mesh, const vector<double>& u0,
                                                    const vector<int>& target, const Config& cfg) {
    DegreeDiagnostics dd;
    dd.ran = true;
    dd.target_degree = target.empty() ? vector<int>(mesh.N, 6) : target;
    vector<int> degree = vertex_degrees(mesh);
    dd.target_degree_radius_corr = pearson_corr_degree_u(dd.target_degree, u0);
    dd.degree_radius_corr = pearson_corr_degree_u(degree, u0);
    double sumsq = 0.0;
    int matches = 0;
    dd.max_degree_error = 0;
    for (int i = 0; i < mesh.N; ++i) {
        int err = degree[i] - dd.target_degree[i];
        sumsq += static_cast<double>(err * err);
        dd.max_degree_error = max(dd.max_degree_error, abs(err));
        if (err == 0) matches++;
        dd.target_hist[dd.target_degree[i]]++;
    }
    dd.E_deg = mesh.N > 0 ? sumsq / mesh.N : 0.0;
    dd.degree_match_fraction = mesh.N > 0 ? static_cast<double>(matches) / mesh.N : 0.0;
    double corr_penalty = max(0.0, cfg.corr_target - dd.degree_radius_corr);
    dd.S_degree = cfg.w_deg * dd.E_deg +
                  cfg.w_degmax * static_cast<double>(dd.max_degree_error) +
                  cfg.w_corr * corr_penalty * corr_penalty;
    return dd;
}

static double linear_total_score(const RunOutput& out, const Config& cfg) {
    double S = out.linear.ran ? out.linear.S_pred : 1e100;
    if (cfg.linear_surgery_use_degree && out.degree.ran && isfinite(out.degree.S_degree)) {
        S += cfg.alpha_degree * out.degree.S_degree;
    }
    return S;
}

static double predictor_score(const LinearResponseResult& lr, const Config& cfg) {
    auto logterm = [](double x, double scale) {
        scale = max(scale, 1e-300);
        return log10(1.0 + max(0.0, x) / scale);
    };
    double corr = isfinite(lr.degree_radius_corr) ? lr.degree_radius_corr : 0.0;
    double degree_radius_mismatch = max(0.0, cfg.corr_target - corr);
    degree_radius_mismatch *= degree_radius_mismatch;
    return cfg.w_res * logterm(lr.rms_K_lin, cfg.K_rms_scale) +
           cfg.w_max * logterm(lr.max_K_lin, cfg.K_max_scale) +
           cfg.w_du * pow(lr.du_rms / max(cfg.du_rms_allowed, 1e-300), 2.0) +
           cfg.w_dumax * pow(lr.du_max / max(cfg.du_max_allowed, 1e-300), 2.0) +
           cfg.w_sing * (lr.singular_response_flag ? 1.0 : 0.0) +
           cfg.w_mis * degree_radius_mismatch;
}

static LinearResponseResult compute_linear_response(const Mesh& mesh, const vector<double>& u0, const Config& cfg) {
    LinearResponseResult lr;
    lr.ran = true;
    lr.status = mesh.N > 1024 ? "warning: dense finite-difference L is expensive for N > 1024" : "ok";
    const int N = mesh.N;
    if (N == 0 || static_cast<int>(u0.size()) != N) {
        lr.status = "invalid u0 size";
        return lr;
    }

    State base;
    compute_state_inplace(mesh, u0, base);
    lr.K0 = base.K;
    lr.rms_K0 = base.rms_K;
    lr.max_K0 = base.max_abs_K;

    vector<double> L(size_t(N) * N, 0.0);
    vector<double> up = u0, um = u0;
    State sp, sm;
    double h = max(cfg.lin_h, 1e-12);
    for (int j = 0; j < N; ++j) {
        up[j] += h;
        um[j] -= h;
        compute_state_inplace(mesh, up, sp);
        compute_state_inplace(mesh, um, sm);
        for (int i = 0; i < N; ++i) {
            L[size_t(i) * N + j] = (sp.K[i] - sm.K[i]) / (2.0 * h);
        }
        up[j] = u0[j];
        um[j] = u0[j];
    }

    vector<double> rhs(N, 0.0);
    for (int j = 0; j < N; ++j) {
        double s = 0.0;
        for (int i = 0; i < N; ++i) s += L[size_t(i) * N + j] * lr.K0[i];
        rhs[j] = -s;
    }

    auto apply_A = [&](const vector<double>& x, vector<double>& y) {
        vector<double> tmp(N, 0.0);
        for (int i = 0; i < N; ++i) {
            double s = 0.0;
            const size_t row = size_t(i) * N;
            for (int j = 0; j < N; ++j) s += L[row + j] * x[j];
            tmp[i] = s;
        }
        y.assign(N, 0.0);
        for (int j = 0; j < N; ++j) {
            double s = cfg.lin_lambda * x[j];
            for (int i = 0; i < N; ++i) s += L[size_t(i) * N + j] * tmp[i];
            y[j] = s;
        }
    };

    lr.du.assign(N, 0.0);
    vector<double> r = rhs;
    vector<double> p = r;
    vector<double> Ap(N, 0.0);
    double rsold = inner_product(r.begin(), r.end(), r.begin(), 0.0);
    double rhs_norm = sqrt(max(rsold, 0.0));
    if (rhs_norm > 0.0) {
        double tol2 = pow(max(cfg.lin_cg_tol, 1e-16) * rhs_norm, 2.0);
        for (int it = 1; it <= cfg.lin_cg_max_iter; ++it) {
            apply_A(p, Ap);
            double denom = inner_product(p.begin(), p.end(), Ap.begin(), 0.0);
            if (fabs(denom) < 1e-300) break;
            double alpha = rsold / denom;
            for (int i = 0; i < N; ++i) {
                lr.du[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }
            double rsnew = inner_product(r.begin(), r.end(), r.begin(), 0.0);
            lr.cg_iters = it;
            if (rsnew <= tol2) {
                rsold = rsnew;
                break;
            }
            double beta = rsnew / max(rsold, 1e-300);
            for (int i = 0; i < N; ++i) p[i] = r[i] + beta * p[i];
            rsold = rsnew;
        }
    }
    lr.cg_residual = sqrt(max(rsold, 0.0));

    lr.K_lin.assign(N, 0.0);
    lr.rms_K_lin = 0.0;
    lr.max_K_lin = 0.0;
    for (int i = 0; i < N; ++i) {
        double s = lr.K0[i];
        const size_t row = size_t(i) * N;
        for (int j = 0; j < N; ++j) s += L[row + j] * lr.du[j];
        lr.K_lin[i] = s;
        lr.rms_K_lin += s * s;
        lr.max_K_lin = max(lr.max_K_lin, fabs(s));
    }
    lr.rms_K_lin = sqrt(lr.rms_K_lin / static_cast<double>(N));

    lr.du_l2 = 0.0;
    lr.du_max = 0.0;
    for (double x : lr.du) {
        lr.du_l2 += x * x;
        lr.du_max = max(lr.du_max, fabs(x));
    }
    lr.du_l2 = sqrt(lr.du_l2);
    lr.du_rms = lr.du_l2 / sqrt(static_cast<double>(N));
    lr.response_ratio = lr.rms_K_lin / max(lr.rms_K0, 1e-300);
    lr.singular_response_flag = lr.du_max > cfg.du_max_allowed || lr.du_rms > cfg.du_rms_allowed;

    vector<int> degree = vertex_degrees(mesh);
    lr.degree_radius_corr = pearson_corr_degree_u(degree, u0);

    double mean_u = accumulate(u0.begin(), u0.end(), 0.0) / N;
    double mean_d = accumulate(degree.begin(), degree.end(), 0.0) / static_cast<double>(N);
    double cov = 0.0, varu = 0.0;
    for (int i = 0; i < N; ++i) {
        double du = u0[i] - mean_u;
        cov += du * (degree[i] - mean_d);
        varu += du * du;
    }
    double beta = varu > 1e-30 ? cov / varu : 0.0;
    double alpha = mean_d - beta * mean_u;
    lr.vertex_mismatch.assign(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double expected_d = alpha + beta * u0[i];
        double local_degree_mismatch = fabs(static_cast<double>(degree[i]) - expected_d);
        lr.vertex_mismatch[i] = fabs(lr.K_lin[i]) +
                                cfg.mismatch_a * fabs(lr.du[i]) +
                                cfg.mismatch_b * local_degree_mismatch;
    }
    lr.edge_mismatch.assign(mesh.edges.size(), 0.0);
    for (int ei = 0; ei < static_cast<int>(mesh.edges.size()); ++ei) {
        const EdgeAdj& e = mesh.edges[ei];
        lr.edge_mismatch[ei] = lr.vertex_mismatch[e.a] + lr.vertex_mismatch[e.b] +
                               cfg.mismatch_c * fabs(lr.du[e.a] - lr.du[e.b]);
    }

    lr.S_pred = predictor_score(lr, cfg);
    lr.linear_predicted_good = lr.rms_K_lin < cfg.predictor_rms_tol &&
                               lr.du_rms < cfg.du_rms_allowed &&
                               lr.du_max < cfg.du_max_allowed;
    return lr;
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
                                                      const DevelopResult& dr, const Config& cfg,
                                                      const DistributionDiagnostics& dd) {
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
        pd.fixed_intrinsic_good = dd.fixed_intrinsic_good;
        pd.weak_distribution_preserved = dd.weak_distribution_preserved;
        pd.distribution_preserved = dd.distribution_preserved;
        pd.physical_INC_like_candidate_v6 = false;
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
        bool bad_len = false;
        for (int e = 0; e < 3; ++e) {
            int q = (e + 1) % 3;
            double len = sqrt(dist2(df.x[e], df.y[e], df.x[q], df.y[q]));
            if (!isfinite(len)) bad_len = true;
        }
        if (fabs(area) < cfg.area_tol || bad_len || !isfinite(area)) pd.face_degenerate_count++;
        if (area < -cfg.area_tol) pd.face_negative_signed_count++;
        pd.face_min_area = min(pd.face_min_area, area);
        pd.face_max_area = max(pd.face_max_area, area);
        area_sum += area;
        area_count++;
    }
    pd.face_bad_area_count = pd.face_degenerate_count;
    if (area_count > 0) {
        pd.face_mean_area = area_sum / area_count;
    } else {
        pd.face_min_area = pd.face_max_area = pd.face_mean_area = numeric_limits<double>::quiet_NaN();
    }

    if (cfg.check_orientation_neighbors) {
        for (const EdgeAdj& edge : mesh.edges) {
            if (edge.f0 < 0 || edge.f1 < 0 || !dr.faces[edge.f0].done || !dr.faces[edge.f1].done) continue;
            const Face& f0 = mesh.faces[edge.f0];
            const Face& f1 = mesh.faces[edge.f1];
            const DevelopedFace& d0 = dr.faces[edge.f0];
            const DevelopedFace& d1 = dr.faces[edge.f1];
            int p0 = edge.e0;
            int q0 = (edge.e0 + 1) % 3;
            int o0 = (edge.e0 + 2) % 3;
            int o1 = opposite_index(f1, edge.a, edge.b);
            int a0 = corner_index_for_vertex(f0, edge.a);
            int b0 = corner_index_for_vertex(f0, edge.b);
            int a1 = corner_index_for_vertex(f1, edge.a);
            int b1 = corner_index_for_vertex(f1, edge.b);
            if (o1 < 0 || a0 < 0 || b0 < 0 || a1 < 0 || b1 < 0) continue;
            double offx = 0.5 * ((d0.x[a0] - d1.x[a1]) + (d0.x[b0] - d1.x[b1]));
            double offy = 0.5 * ((d0.y[a0] - d1.y[a1]) + (d0.y[b0] - d1.y[b1]));
            double sx0 = orient2d(d0.x[p0], d0.y[p0], d0.x[q0], d0.y[q0], d0.x[o0], d0.y[o0]);
            double sx1 = orient2d(d0.x[p0], d0.y[p0], d0.x[q0], d0.y[q0], d1.x[o1] + offx, d1.y[o1] + offy);
            if (fabs(sx0) > cfg.area_tol && fabs(sx1) > cfg.area_tol && sx0 * sx1 > 0.0) {
                pd.orientation_neighbor_violation_count++;
            }
        }
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

    double edge_mean = 0.0;
    double edge_s2 = 0.0;
    for (const EdgeError& e : pd.edge_errors) {
        edge_mean += e.dist;
        edge_s2 += e.dist * e.dist;
    }
    if (!pd.edge_errors.empty()) {
        edge_mean /= static_cast<double>(pd.edge_errors.size());
        double var = max(0.0, edge_s2 / static_cast<double>(pd.edge_errors.size()) - edge_mean * edge_mean);
        pd.edge_length_cv = edge_mean > 0.0 ? sqrt(var) / edge_mean : numeric_limits<double>::infinity();
    }
    vector<int> degree(mesh.N, 0);
    for (const EdgeAdj& e : mesh.edges) {
        degree[e.a]++;
        degree[e.b]++;
    }
    int degree6 = 0;
    for (int d : degree) {
        if (d == 6) degree6++;
    }
    pd.degree6_fraction = mesh.N > 0 ? static_cast<double>(degree6) / mesh.N : numeric_limits<double>::quiet_NaN();
    pd.noncrystalline_topology_heuristic = pd.degree6_fraction < 0.95;

    pd.topology_geometry_valid = pd.face_degenerate_count == 0 &&
                                 (!cfg.check_crossings || pd.edge_crossing_count == 0) &&
                                 (!cfg.require_consistent_orientation || pd.face_negative_signed_count == 0) &&
                                 (!cfg.check_orientation_neighbors || pd.orientation_neighbor_violation_count == 0);
    pd.valid_disk_packing_candidate = pd.intrinsic_endpoint &&
                                      pd.developed_consistent &&
                                      pd.contact_edges_valid &&
                                      pd.nonedge_overlap_free &&
                                      pd.topology_geometry_valid;
    pd.strong_INC_candidate = pd.valid_disk_packing_candidate && pd.local_theta_valid;
    pd.physical_INC_like_candidate = pd.strong_INC_candidate && pd.physically_reasonable_radius;
    pd.fixed_intrinsic_good = dd.fixed_intrinsic_good;
    pd.weak_distribution_preserved = dd.weak_distribution_preserved;
    pd.distribution_preserved = dd.distribution_preserved;
    pd.physical_INC_like_candidate_v6 = pd.strong_INC_candidate &&
                                       pd.distribution_preserved &&
                                       pd.physically_reasonable_radius;
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

static vector<string> parse_string_list(const string& s) {
    vector<string> out;
    string item;
    stringstream ss(s);
    while (getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
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
        else if (arg == "--require_consistent_orientation") cfg.require_consistent_orientation = parse_bool_int(val(arg));
        else if (arg == "--check_orientation_neighbors") cfg.check_orientation_neighbors = parse_bool_int(val(arg));
        else if (arg == "--check_crossings") {
            cfg.check_crossings = parse_bool_int(val(arg));
            cfg.check_crossings_explicit = true;
        }
        else if (arg == "--intersect_tol") cfg.intersect_tol = stod(val(arg));
        else if (arg == "--hist_bins") cfg.hist_bins = stoi(val(arg));
        else if (arg == "--target_max_ratio") cfg.target_max_ratio = stod(val(arg));
        else if (arg == "--target_delta") cfg.target_delta = stod(val(arg));
        else if (arg == "--radius_mode") cfg.radius_mode = val(arg);
        else if (arg == "--target_dist") cfg.target_dist = val(arg);
        else if (arg == "--target_radius_ratio") cfg.target_radius_ratio = stod(val(arg));
        else if (arg == "--radius_seed") cfg.radius_seed = stoull(val(arg));
        else if (arg == "--radius_noise_std") cfg.radius_noise_std = stod(val(arg));
        else if (arg == "--fixed_tol_K") cfg.fixed_tol_K = stod(val(arg));
        else if (arg == "--fixed_tol_rms_K") cfg.fixed_tol_rms_K = stod(val(arg));
        else if (arg == "--lambda_u") cfg.lambda_u = stod(val(arg));
        else if (arg == "--lambda_sort") cfg.lambda_sort = stod(val(arg));
        else if (arg == "--lambda_quantile") cfg.lambda_quantile = stod(val(arg));
        else if (arg == "--lambda_delta") cfg.lambda_delta = stod(val(arg));
        else if (arg == "--lambda_ratio") cfg.lambda_ratio = stod(val(arg));
        else if (arg == "--dist_tol_sort") cfg.dist_tol_sort = stod(val(arg));
        else if (arg == "--dist_tol_delta") cfg.dist_tol_delta = stod(val(arg));
        else if (arg == "--dist_tol_ratio") cfg.dist_tol_ratio = stod(val(arg));
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
        else if (arg == "--surgery_batch_flips") cfg.surgery_batch_flips = stoi(val(arg));
        else if (arg == "--proposal_mode") cfg.proposal_mode = val(arg);
        else if (arg == "--K_scale") cfg.K_scale = stod(val(arg));
        else if (arg == "--rmsK_scale") cfg.rmsK_scale = stod(val(arg));
        else if (arg == "--dev_scale") cfg.dev_scale = stod(val(arg));
        else if (arg == "--period_scale") cfg.period_scale = stod(val(arg));
        else if (arg == "--overlap_scale") cfg.overlap_scale = stod(val(arg));
        else if (arg == "--wK") cfg.wK = stod(val(arg));
        else if (arg == "--wRmsK") cfg.wRmsK = stod(val(arg));
        else if (arg == "--wDev") cfg.wDev = stod(val(arg));
        else if (arg == "--wPer") cfg.wPer = stod(val(arg));
        else if (arg == "--wOv") cfg.wOv = stod(val(arg));
        else if (arg == "--wOvN") cfg.wOvN = stod(val(arg));
        else if (arg == "--wCross") cfg.wCross = stod(val(arg));
        else if (arg == "--wDeg") cfg.wDeg = stod(val(arg));
        else if (arg == "--wOrient") cfg.wOrient = stod(val(arg));
        else if (arg == "--wRad") cfg.wRad = stod(val(arg));
        else if (arg == "--wDelta") cfg.wDelta = stod(val(arg));
        else if (arg == "--wSort") cfg.wSort = stod(val(arg));
        else if (arg == "--wRatio") cfg.wRatio = stod(val(arg));
        else if (arg == "--multistart_surgery") cfg.multistart_surgery = parse_bool_int(val(arg));
        else if (arg == "--multistart_count") cfg.multistart_count = stoi(val(arg));
        else if (arg == "--multistart_out") cfg.multistart_out = val(arg);
        else if (arg == "--pareto_scan") cfg.pareto_scan = parse_bool_int(val(arg));
        else if (arg == "--pareto_out") cfg.pareto_out = val(arg);
        else if (arg == "--pareto_delta_list") cfg.pareto_delta_list = parse_double_list(val(arg));
        else if (arg == "--pareto_ratio_list") cfg.pareto_ratio_list = parse_double_list(val(arg));
        else if (arg == "--pareto_modes") cfg.pareto_modes = parse_string_list(val(arg));
        else if (arg == "--pareto_seeds") cfg.pareto_seeds = parse_seed_list(val(arg));
        else if (arg == "--lin_h") cfg.lin_h = stod(val(arg));
        else if (arg == "--lin_lambda") cfg.lin_lambda = stod(val(arg));
        else if (arg == "--lin_cg_tol") cfg.lin_cg_tol = stod(val(arg));
        else if (arg == "--lin_cg_max_iter") cfg.lin_cg_max_iter = stoi(val(arg));
        else if (arg == "--du_max_allowed") cfg.du_max_allowed = stod(val(arg));
        else if (arg == "--du_rms_allowed") cfg.du_rms_allowed = stod(val(arg));
        else if (arg == "--K_rms_scale") cfg.K_rms_scale = stod(val(arg));
        else if (arg == "--K_max_scale") cfg.K_max_scale = stod(val(arg));
        else if (arg == "--w_res") cfg.w_res = stod(val(arg));
        else if (arg == "--w_max") cfg.w_max = stod(val(arg));
        else if (arg == "--w_du") cfg.w_du = stod(val(arg));
        else if (arg == "--w_dumax") cfg.w_dumax = stod(val(arg));
        else if (arg == "--w_sing") cfg.w_sing = stod(val(arg));
        else if (arg == "--w_mis") cfg.w_mis = stod(val(arg));
        else if (arg == "--corr_target") cfg.corr_target = stod(val(arg));
        else if (arg == "--mismatch_a") cfg.mismatch_a = stod(val(arg));
        else if (arg == "--mismatch_b") cfg.mismatch_b = stod(val(arg));
        else if (arg == "--mismatch_c") cfg.mismatch_c = stod(val(arg));
        else if (arg == "--linear_surgery_steps") cfg.linear_surgery_steps = stoi(val(arg));
        else if (arg == "--linear_surgery_trials_per_step") cfg.linear_surgery_trials_per_step = stoi(val(arg));
        else if (arg == "--linear_surgery_temperature") cfg.linear_surgery_temperature = stod(val(arg));
        else if (arg == "--refine_after_linear") cfg.refine_after_linear = parse_bool_int(val(arg));
        else if (arg == "--refine_max_iter") cfg.refine_max_iter = stoi(val(arg));
        else if (arg == "--refine_tol") cfg.refine_tol = stod(val(arg));
        else if (arg == "--refine_du_bound") cfg.refine_du_bound = stod(val(arg));
        else if (arg == "--predictor_rms_tol") cfg.predictor_rms_tol = stod(val(arg));
        else if (arg == "--screen_topologies") cfg.screen_topologies = parse_bool_int(val(arg));
        else if (arg == "--screen_count") cfg.screen_count = stoi(val(arg));
        else if (arg == "--screen_flips") cfg.screen_flips = stoi(val(arg));
        else if (arg == "--screen_out") cfg.screen_out = val(arg);
        else if (arg == "--screen_same_radii") cfg.screen_same_radii = parse_bool_int(val(arg));
        else if (arg == "--linear_pareto_scan") cfg.linear_pareto_scan = parse_bool_int(val(arg));
        else if (arg == "--degree_target_mode") cfg.degree_target_mode = val(arg);
        else if (arg == "--frac_deg4") cfg.frac_deg4 = stod(val(arg));
        else if (arg == "--frac_deg5") cfg.frac_deg5 = stod(val(arg));
        else if (arg == "--frac_deg6") cfg.frac_deg6 = stod(val(arg));
        else if (arg == "--frac_deg7") cfg.frac_deg7 = stod(val(arg));
        else if (arg == "--frac_deg8") cfg.frac_deg8 = stod(val(arg));
        else if (arg == "--w_deg") cfg.w_deg = stod(val(arg));
        else if (arg == "--w_degmax") cfg.w_degmax = stod(val(arg));
        else if (arg == "--w_corr") cfg.w_corr = stod(val(arg));
        else if (arg == "--degree_precondition") cfg.degree_precondition = parse_bool_int(val(arg));
        else if (arg == "--degree_pre_steps") cfg.degree_pre_steps = stoi(val(arg));
        else if (arg == "--degree_pre_trials") cfg.degree_pre_trials = stoi(val(arg));
        else if (arg == "--degree_pre_temperature") cfg.degree_pre_temperature = stod(val(arg));
        else if (arg == "--degree_proposal_mode") cfg.degree_proposal_mode = val(arg);
        else if (arg == "--alpha_degree") cfg.alpha_degree = stod(val(arg));
        else if (arg == "--linear_surgery_use_degree") cfg.linear_surgery_use_degree = parse_bool_int(val(arg));
        else if (arg == "--beta_degree") cfg.beta_degree = stod(val(arg));
        else if (arg == "--degree_fraction_scan") cfg.degree_fraction_scan = parse_bool_int(val(arg));
        else if (arg == "--degree_fraction_out") cfg.degree_fraction_out = val(arg);
        else if (arg == "--compare_topology_strategies") cfg.compare_topology_strategies = parse_bool_int(val(arg));
        else if (arg == "--strategy_out") cfg.strategy_out = val(arg);
        else if (arg == "--test") cfg.test = true;
        else if (arg == "-h" || arg == "--help") {
            cout << "Usage: ./inc_ricci_v8 [--nx N] [--ny N] [--flips N] [--seed N] [--degree_target_mode quantile]\n";
            exit(0);
        } else {
            throw runtime_error("unknown argument: " + arg);
        }
    }
    if (cfg.nx <= 2 || cfg.ny <= 2) throw runtime_error("nx and ny must be > 2");
    if (cfg.flips < 0) throw runtime_error("flips must be >= 0");
    if (cfg.max_iter < 0) throw runtime_error("max_iter must be >= 0");
    if (cfg.tol <= 0.0) throw runtime_error("tol must be positive");
    if (cfg.method != "ricci") throw runtime_error("v8 production path supports --method ricci");
    if (cfg.u_min >= cfg.u_max) throw runtime_error("u_min must be < u_max");
    if (cfg.validate_every < 0) throw runtime_error("validate_every must be >= 0");
    if (cfg.progress_every < 0) throw runtime_error("progress_every must be >= 0");
    if (cfg.threads <= 0) cfg.threads = 1;
    if (cfg.hist_bins <= 0) throw runtime_error("hist_bins must be > 0");
    if (cfg.overlap_mode != "full" && cfg.overlap_mode != "cell") throw runtime_error("overlap_mode must be full or cell");
    if (cfg.radius_mode != "free" && cfg.radius_mode != "fixed" && cfg.radius_mode != "weak") {
        throw runtime_error("radius_mode must be free, fixed, or weak");
    }
    if (cfg.target_dist != "lognormal" && cfg.target_dist != "uniform_u" && cfg.target_dist != "powerlaw_like") {
        throw runtime_error("target_dist must be lognormal, uniform_u, or powerlaw_like");
    }
    if (cfg.target_radius_ratio <= 1.0) throw runtime_error("target_radius_ratio must be > 1");
    if (cfg.surgery_metric != "overlap" && cfg.surgery_metric != "radius" && cfg.surgery_metric != "mixed" &&
        cfg.surgery_metric != "v5" && cfg.surgery_metric != "v6" && cfg.surgery_metric != "linear") {
        throw runtime_error("surgery_metric must be overlap, radius, mixed, v5, v6, or linear");
    }
    if (cfg.surgery_batch_flips <= 0) throw runtime_error("surgery_batch_flips must be > 0");
    if (cfg.proposal_mode != "random" && cfg.proposal_mode != "highK" &&
        cfg.proposal_mode != "overlap" && cfg.proposal_mode != "mismatch") {
        throw runtime_error("proposal_mode must be random, highK, overlap, or mismatch");
    }
    if (cfg.multistart_count <= 0) throw runtime_error("multistart_count must be > 0");
    if (cfg.lin_h <= 0.0) throw runtime_error("lin_h must be positive");
    if (cfg.lin_lambda <= 0.0) throw runtime_error("lin_lambda must be positive");
    if (cfg.linear_surgery_steps < 0) throw runtime_error("linear_surgery_steps must be >= 0");
    if (cfg.linear_surgery_trials_per_step <= 0) throw runtime_error("linear_surgery_trials_per_step must be > 0");
    if (cfg.screen_count <= 0) throw runtime_error("screen_count must be > 0");
    if (cfg.screen_flips < 0) throw runtime_error("screen_flips must be >= 0");
    if (cfg.degree_target_mode != "quantile" && cfg.degree_target_mode != "linear" && cfg.degree_target_mode != "none") {
        throw runtime_error("degree_target_mode must be quantile, linear, or none");
    }
    if (cfg.degree_proposal_mode != "random" && cfg.degree_proposal_mode != "greedy" &&
        cfg.degree_proposal_mode != "annealed") {
        throw runtime_error("degree_proposal_mode must be random, greedy, or annealed");
    }
    if (cfg.degree_pre_steps < 0 || cfg.degree_pre_trials <= 0) throw runtime_error("invalid degree precondition settings");
    if (cfg.pareto_delta_list.empty()) cfg.pareto_delta_list = {0.20, 0.25, 0.30, 0.35, 0.40};
    if (cfg.pareto_ratio_list.empty()) cfg.pareto_ratio_list = {3.0, 5.0, 8.0, 12.0};
    if (cfg.pareto_modes.empty()) cfg.pareto_modes = {"fixed", "weak"};
    if (cfg.pareto_seeds.empty()) cfg.pareto_seeds = {1, 2, 3};
    if (cfg.surgery_metric == "linear" && cfg.radius_mode == "free") cfg.radius_mode = "fixed";
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

static bool attempt_indexed_flip_with_info(Mesh& mesh, int edge_idx, int& a, int& b, int& c, int& d) {
    unordered_map<uint64_t, int> edge_index;
    string error;
    if (!build_adjacency(mesh, &error, &edge_index)) return false;
    if (edge_idx < 0 || edge_idx >= static_cast<int>(mesh.edges.size())) return false;
    const EdgeAdj edge = mesh.edges[edge_idx];
    const Face& f0 = mesh.faces[edge.f0];
    const Face& f1 = mesh.faces[edge.f1];
    int o0 = opposite_index(f0, edge.a, edge.b);
    int o1 = opposite_index(f1, edge.a, edge.b);
    if (o0 < 0 || o1 < 0) return false;
    a = edge.a;
    b = edge.b;
    c = f0.c[o0].v;
    d = f1.c[o1].v;
    return attempt_flip_edge_index(mesh, edge_idx, edge_index);
}

static int run_degree_precondition(Mesh& mesh, const vector<double>& u0, const vector<int>& target,
                                   const Config& cfg, const filesystem::path& out_dir) {
    filesystem::create_directories(out_dir);
    ofstream stats(out_dir / "degree_precondition_stats.csv");
    stats << "step,accepted,S_degree,E_deg,max_degree_error,degree_match_fraction,degree_radius_corr,"
             "flipped_a,flipped_b,flipped_c,flipped_d\n";
    stats << scientific << setprecision(17);

    mt19937_64 rng(cfg.seed + 0xd36d5eedULL);
    uniform_real_distribution<double> uni(0.0, 1.0);
    DegreeDiagnostics current = compute_degree_diagnostics(mesh, u0, target, cfg);
    int accepts = 0;

    for (int step = 1; step <= cfg.degree_pre_steps; ++step) {
        Mesh best_mesh;
        DegreeDiagnostics best_diag;
        double best_score = numeric_limits<double>::infinity();
        bool have_best = false;
        int ba = -1, bb = -1, bc = -1, bd = -1;

        string error;
        if (!build_adjacency(mesh, &error)) break;
        uniform_int_distribution<int> pick(0, max(0, static_cast<int>(mesh.edges.size()) - 1));

        int trials = cfg.degree_proposal_mode == "random" ? 1 : cfg.degree_pre_trials;
        for (int t = 0; t < trials; ++t) {
            Mesh proposal = mesh;
            if (!build_adjacency(proposal, &error)) continue;
            int edge_idx = pick(rng);
            int a = -1, b = -1, c = -1, d = -1;
            if (!attempt_indexed_flip_with_info(proposal, edge_idx, a, b, c, d)) continue;
            if (!validate_mesh(proposal, &error)) continue;
            DegreeDiagnostics diag = compute_degree_diagnostics(proposal, u0, target, cfg);
            if (diag.S_degree < best_score) {
                best_score = diag.S_degree;
                best_mesh = std::move(proposal);
                best_diag = std::move(diag);
                ba = a; bb = b; bc = c; bd = d;
                have_best = true;
            }
        }

        bool accepted = false;
        if (have_best && best_score < current.S_degree) {
            accepted = true;
        } else if (have_best && (cfg.degree_proposal_mode == "annealed" || cfg.degree_pre_temperature > 0.0)) {
            double T = max(cfg.degree_pre_temperature, 1e-12);
            double p = exp(-(best_score - current.S_degree) / T);
            accepted = uni(rng) < p;
        }
        if (accepted) {
            mesh = std::move(best_mesh);
            current = std::move(best_diag);
            accepts++;
        } else {
            ba = bb = bc = bd = -1;
        }

        stats << step << ',' << (accepted ? 1 : 0) << ',' << current.S_degree << ','
              << current.E_deg << ',' << current.max_degree_error << ','
              << current.degree_match_fraction << ',' << current.degree_radius_corr << ','
              << ba << ',' << bb << ',' << bc << ',' << bd << '\n';
    }
    return accepts;
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

        vector<double> u0;
        vector<double> u;
        if (cfg.radius_mode == "free") {
            u = random_log_radii(out.mesh.N, rng, 0.1);
        } else {
            u0 = generate_target_log_radii(out.mesh.N, cfg);
            u = u0;
        }
        vector<int> target_degree;
        if (!u0.empty()) {
            target_degree = assign_target_degrees(u0, cfg);
            if (cfg.degree_precondition) {
                int pre_accepts = run_degree_precondition(out.mesh, u0, target_degree, cfg, filesystem::path(cfg.out));
                out.flips_accepted += pre_accepts;
                if (!validate_mesh(out.mesh, &error)) throw runtime_error("mesh invalid after degree precondition: " + error);
            }
            out.degree = compute_degree_diagnostics(out.mesh, u0, target_degree, cfg);
        }
        out.time.topology = now_sec() - t0;
        t0 = now_sec();
        out.ricci = run_ricci(out.mesh, std::move(u), cfg, u0, write_files, verbose);
        out.time.ricci = now_sec() - t0;
        out.dist = compute_distribution_diagnostics(out.ricci.state.u, u0, cfg,
                                                    out.ricci.state.max_abs_K,
                                                    out.ricci.state.rms_K);
        if (!u0.empty()) {
            vector<int> target = assign_target_degrees(u0, cfg);
            out.degree = compute_degree_diagnostics(out.mesh, u0, target, cfg);
        }
        if ((cfg.surgery_metric == "linear" || cfg.linear_pareto_scan || cfg.screen_topologies) && !out.dist.u0.empty()) {
            out.linear = compute_linear_response(out.mesh, out.dist.u0, cfg);
        }

        if (cfg.develop) {
            t0 = now_sec();
            out.develop = develop_mesh(out.mesh, out.ricci.state.r);
            out.time.develop = now_sec() - t0;
            out.packing = compute_packing_diagnostics(out.mesh, out.ricci, out.develop, cfg, out.dist);
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
    const DistributionDiagnostics& dd = out.dist;
    const DegreeDiagnostics& gd = out.degree;

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
    if (!dd.u0.empty()) {
        ofstream f(filesystem::path(cfg.out) / "target_radii.csv");
        f << "vertex,u0,r0\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < mesh.N; ++i) f << i << ',' << dd.u0[i] << ',' << dd.r0[i] << '\n';
    }
    if (gd.ran && !dd.u0.empty()) {
        ofstream f(filesystem::path(cfg.out) / "target_degrees.csv");
        f << "vertex,u0,r0,d_target\n";
        f << scientific << setprecision(17);
        for (int i = 0; i < mesh.N; ++i) {
            f << i << ',' << dd.u0[i] << ',' << dd.r0[i] << ',' << gd.target_degree[i] << '\n';
        }
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
            f << "fixed_intrinsic_good: " << (pd.fixed_intrinsic_good ? "yes" : "no") << "\n";
            f << "weak_distribution_preserved: " << (pd.weak_distribution_preserved ? "yes" : "no") << "\n";
            f << "distribution_preserved: " << (pd.distribution_preserved ? "yes" : "no") << "\n";
            f << "physical_INC_like_candidate_v6: " << (pd.physical_INC_like_candidate_v6 ? "yes" : "no") << "\n";
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
            f << "face_degenerate_count: " << pd.face_degenerate_count << "\n";
            f << "face_negative_signed_count: " << pd.face_negative_signed_count << "\n";
            f << "orientation_neighbor_violation_count: " << pd.orientation_neighbor_violation_count << "\n";
            f << "edge_crossing_count: " << pd.edge_crossing_count << "\n";
            f << "edge_crossing_fraction: " << pd.edge_crossing_fraction << "\n";
            f << "cell_skew_sin: " << pd.cell_skew_sin << "\n";
            f << "cell_skew_warning: " << (pd.cell_skew_warning ? "yes" : "no") << "\n";
        }
        {
            ofstream f(filesystem::path(cfg.out) / "disorder_stats.txt");
            f << scientific << setprecision(17);
            f << "degree6_fraction: " << pd.degree6_fraction << "\n";
            f << "degree_histogram:";
            for (const auto& kv : degree_histogram(mesh)) f << " " << kv.first << ":" << kv.second;
            f << "\n";
            f << "edge_length_cv: " << pd.edge_length_cv << "\n";
            f << "radius_cv: " << pd.radius_stats.cv_r << "\n";
            f << "noncrystalline_topology_heuristic: " << (pd.noncrystalline_topology_heuristic ? "yes" : "no") << "\n";
        }
        {
            ofstream f(filesystem::path(cfg.out) / "distribution_stats.txt");
            f << scientific << setprecision(17);
            f << "radius_mode: " << cfg.radius_mode << "\n";
            f << "target_dist: " << cfg.target_dist << "\n";
            f << "target_delta: " << cfg.target_delta << "\n";
            f << "target_radius_ratio: " << cfg.target_radius_ratio << "\n";
            f << "D_u: " << dd.D_u << "\n";
            f << "D_sort: " << dd.D_sort << "\n";
            f << "D_quantile: " << dd.D_quantile << "\n";
            f << "D_delta: " << dd.D_delta << "\n";
            f << "D_ratio: " << dd.D_ratio << "\n";
            f << "current_delta: " << dd.current_delta << "\n";
            f << "current_radius_ratio: " << dd.current_radius_ratio << "\n";
            f << "q05_u: " << dd.q05_u << "\n";
            f << "q50_u: " << dd.q50_u << "\n";
            f << "q95_u: " << dd.q95_u << "\n";
            f << "q05_u0: " << dd.q05_u0 << "\n";
            f << "q50_u0: " << dd.q50_u0 << "\n";
            f << "q95_u0: " << dd.q95_u0 << "\n";
            f << "fixed_intrinsic_good: " << (dd.fixed_intrinsic_good ? "yes" : "no") << "\n";
            f << "weak_distribution_preserved: " << (dd.weak_distribution_preserved ? "yes" : "no") << "\n";
            f << "distribution_preserved: " << (dd.distribution_preserved ? "yes" : "no") << "\n";
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
        f << "INC-Ricci v8 degree-radius topology-design run\n";
        f << "N: " << mesh.N << "\n";
        f << "E: " << mesh.edges.size() << "\n";
        f << "F: " << mesh.faces.size() << "\n";
        f << "flips_requested: " << cfg.flips << "\n";
        f << "flips_accepted: " << out.flips_accepted << "\n";
        f << "seed: " << cfg.seed << "\n";
        f << "bounded: " << (cfg.bounded ? "yes" : "no") << "\n";
        f << "radius_mode: " << cfg.radius_mode << "\n";
        f << "target_dist: " << cfg.target_dist << "\n";
        f << "target_delta: " << cfg.target_delta << "\n";
        f << "target_radius_ratio: " << cfg.target_radius_ratio << "\n";
        f << "boundary_rejections: " << rr.boundary_rejections << "\n";
        f << "ricci_iters: " << rr.iterations << "\n";
        f << "final_E_K: " << rr.state.E_K << "\n";
        f << "final_max_abs_K: " << rr.state.max_abs_K << "\n";
        f << "final_rms_K: " << rr.state.rms_K << "\n";
        f << "final_E_reg: " << rr.final_E_reg << "\n";
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
        f << "face_degenerate_count: " << pd.face_degenerate_count << "\n";
        f << "face_negative_signed_count: " << pd.face_negative_signed_count << "\n";
        f << "orientation_neighbor_violation_count: " << pd.orientation_neighbor_violation_count << "\n";
        f << "degree6_fraction: " << pd.degree6_fraction << "\n";
        f << "noncrystalline_topology_heuristic: " << (pd.noncrystalline_topology_heuristic ? "yes" : "no") << "\n";
        f << "valid_disk_packing_candidate: " << (pd.valid_disk_packing_candidate ? "yes" : "no") << "\n";
        f << "strong_INC_candidate: " << (pd.strong_INC_candidate ? "yes" : "no") << "\n";
        f << "physical_INC_like_candidate: " << (pd.physical_INC_like_candidate ? "yes" : "no") << "\n";
        f << "fixed_intrinsic_good: " << (pd.fixed_intrinsic_good ? "yes" : "no") << "\n";
        f << "weak_distribution_preserved: " << (pd.weak_distribution_preserved ? "yes" : "no") << "\n";
        f << "distribution_preserved: " << (pd.distribution_preserved ? "yes" : "no") << "\n";
        f << "physical_INC_like_candidate_v6: " << (pd.physical_INC_like_candidate_v6 ? "yes" : "no") << "\n";
        f << "D_sort: " << dd.D_sort << "\n";
        f << "D_quantile: " << dd.D_quantile << "\n";
        f << "D_delta: " << dd.D_delta << "\n";
        f << "D_ratio: " << dd.D_ratio << "\n";
        f << "S_degree: " << gd.S_degree << "\n";
        f << "E_deg: " << gd.E_deg << "\n";
        f << "max_degree_error: " << gd.max_degree_error << "\n";
        f << "degree_match_fraction: " << gd.degree_match_fraction << "\n";
        f << "degree_radius_corr: " << gd.degree_radius_corr << "\n";
        f << "target_degree_radius_corr: " << gd.target_degree_radius_corr << "\n";
        f << "target_degree_histogram:";
        for (const auto& kv : gd.target_hist) f << " " << kv.first << ":" << kv.second;
        f << "\n";
        f << "time_topology_sec: " << out.time.topology << "\n";
        f << "time_ricci_sec: " << out.time.ricci << "\n";
        f << "time_develop_sec: " << out.time.develop << "\n";
        f << "time_embed_sec: " << out.time.embed << "\n";
        f << "time_total_sec: " << out.time.total << "\n";
        f << "note: v8 tests degree-radius topology design and linear-response prediction; it is not a molecular SMC/FIRE or full packing dynamics solver.\n";
    }
    if (out.linear.ran) {
        const LinearResponseResult& lr = out.linear;
        vector<int> degree = vertex_degrees(mesh);
        {
            ofstream f(filesystem::path(cfg.out) / "linear_response_summary.txt");
            f << scientific << setprecision(17);
            f << "N: " << mesh.N << "\n";
            f << "E: " << mesh.edges.size() << "\n";
            f << "F: " << mesh.faces.size() << "\n";
            f << "target_delta: " << cfg.target_delta << "\n";
            f << "target_radius_ratio: " << cfg.target_radius_ratio << "\n";
            f << "rms_K0: " << lr.rms_K0 << "\n";
            f << "max_K0: " << lr.max_K0 << "\n";
            f << "rms_K_lin: " << lr.rms_K_lin << "\n";
            f << "max_K_lin: " << lr.max_K_lin << "\n";
            f << "du_rms: " << lr.du_rms << "\n";
            f << "du_max: " << lr.du_max << "\n";
            f << "du_l2: " << lr.du_l2 << "\n";
            f << "response_ratio: " << lr.response_ratio << "\n";
            f << "S_pred: " << lr.S_pred << "\n";
            f << "S_degree: " << gd.S_degree << "\n";
            f << "S_total: " << linear_total_score(out, cfg) << "\n";
            f << "degree_radius_corr: " << lr.degree_radius_corr << "\n";
            f << "singular_response_flag: " << (lr.singular_response_flag ? "yes" : "no") << "\n";
            f << "linear_predicted_good: " << (lr.linear_predicted_good ? "yes" : "no") << "\n";
            f << "cg_iters: " << lr.cg_iters << "\n";
            f << "cg_residual: " << lr.cg_residual << "\n";
            f << "degree_histogram:";
            for (const auto& kv : degree_histogram(mesh)) f << " " << kv.first << ":" << kv.second;
            f << "\n";
            map<int, vector<double>> by_degree;
            for (int i = 0; i < mesh.N; ++i) by_degree[degree[i]].push_back(dd.u0.empty() ? rr.state.u[i] : dd.u0[i]);
            f << "mean_u_by_degree:";
            for (const auto& kv : by_degree) {
                double m = accumulate(kv.second.begin(), kv.second.end(), 0.0) / kv.second.size();
                f << " " << kv.first << ":" << m;
            }
            f << "\nvar_u_by_degree:";
            for (const auto& kv : by_degree) {
                double m = accumulate(kv.second.begin(), kv.second.end(), 0.0) / kv.second.size();
                double v = 0.0;
                for (double x : kv.second) v += (x - m) * (x - m);
                v /= kv.second.size();
                f << " " << kv.first << ":" << v;
            }
            f << "\n";
        }
        {
            ofstream f(filesystem::path(cfg.out) / "linear_du.csv");
            f << "vertex,u0,du,u_pred,K0,K_lin,degree,vertex_mismatch\n";
            f << scientific << setprecision(17);
            for (int i = 0; i < mesh.N; ++i) {
                double u0 = dd.u0.empty() ? rr.state.u[i] : dd.u0[i];
                double du = i < static_cast<int>(lr.du.size()) ? lr.du[i] : 0.0;
                f << i << ',' << u0 << ',' << du << ',' << u0 + du << ','
                  << lr.K0[i] << ',' << lr.K_lin[i] << ',' << degree[i] << ','
                  << lr.vertex_mismatch[i] << '\n';
            }
        }
        {
            ofstream f(filesystem::path(cfg.out) / "linear_edges.csv");
            f << "i,j,edge_mismatch\n";
            f << scientific << setprecision(17);
            for (int ei = 0; ei < static_cast<int>(mesh.edges.size()); ++ei) {
                const EdgeAdj& e = mesh.edges[ei];
                f << e.a << ',' << e.b << ',' << lr.edge_mismatch[ei] << '\n';
            }
        }
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
    row.face_degenerate_count = out.packing.face_degenerate_count;
    row.orientation_neighbor_violation_count = out.packing.orientation_neighbor_violation_count;
    row.polydispersity_delta = out.packing.radius_stats.polydispersity_delta;
    row.physically_reasonable_radius = out.packing.physically_reasonable_radius;
    row.valid_disk_packing_candidate = out.packing.valid_disk_packing_candidate;
    row.strong_INC_candidate = out.packing.strong_INC_candidate;
    row.physical_INC_like_candidate = out.packing.physical_INC_like_candidate;
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
         "face_degenerate_count,orientation_neighbor_violation_count,"
         "physically_reasonable_radius,valid_disk_packing_candidate,strong_INC_candidate,physical_INC_like_candidate,"
         "time_total_sec,time_ricci_sec,time_develop_sec,stopped_reason,status\n";
    f << scientific << setprecision(17);
    for (const ScanRow& r : rows) {
        f << r.nx << ',' << r.ny << ',' << r.N << ',' << r.E << ',' << r.F << ','
          << r.flips_requested << ',' << r.flips_accepted << ',' << r.seed << ','
          << r.ricci_iters << ',' << r.final_E_K << ',' << r.final_max_abs_K << ','
          << r.radius_ratio << ',' << r.polydispersity_delta << ',' << r.develop_global_vertex_max_spread << ','
          << r.period_fit_rms << ',' << r.develop_theta_mean_abs << ','
          << r.nonedge_overlap_count << ',' << r.nonedge_max_overlap << ',' << r.edge_crossing_count << ','
          << r.face_bad_area_count << ',' << r.face_degenerate_count << ','
          << r.orientation_neighbor_violation_count << ','
          << (r.physically_reasonable_radius ? 1 : 0) << ','
          << (r.valid_disk_packing_candidate ? 1 : 0) << ','
          << (r.strong_INC_candidate ? 1 : 0) << ','
          << (r.physical_INC_like_candidate ? 1 : 0) << ','
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
    r.face_degenerate_count = out.packing.face_degenerate_count;
    r.orientation_neighbor_violation_count = out.packing.orientation_neighbor_violation_count;
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
         "edge_crossing_count,face_bad_area_count,face_degenerate_count,orientation_neighbor_violation_count,"
         "physically_reasonable_radius,valid_disk_packing_candidate,"
         "strong_INC_candidate,time_total_sec,status\n";
    f << scientific << setprecision(17);
    for (const BoundsRow& r : rows) {
        f << r.nx << ',' << r.ny << ',' << r.seed << ',' << r.flips_requested << ','
          << r.flips_accepted << ',' << r.B << ',' << r.unbounded << ','
          << r.final_E_K << ',' << r.final_max_abs_K << ',' << r.radius_ratio << ','
          << r.polydispersity_delta << ',' << r.boundary_rejections << ',' << r.ricci_iters << ','
          << r.develop_global_vertex_max_spread << ',' << r.period_fit_rms << ','
          << r.nonedge_overlap_count << ',' << r.nonedge_max_overlap << ','
          << r.edge_crossing_count << ',' << r.face_bad_area_count << ',' << r.face_degenerate_count << ','
          << r.orientation_neighbor_violation_count << ','
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

struct ScoreBreakdown {
    double total = 0.0;
    double K = 0.0;
    double rmsK = 0.0;
    double develop = 0.0;
    double period = 0.0;
    double overlap = 0.0;
    double overlap_count = 0.0;
    double crossing = 0.0;
    double degenerate = 0.0;
    double orientation = 0.0;
    double radius = 0.0;
    double delta = 0.0;
    double sort = 0.0;
    double ratio = 0.0;
};

static ScoreBreakdown compute_score_breakdown(const RunOutput& out, const Config& cfg) {
    double maxK = isfinite(out.ricci.state.max_abs_K) ? out.ricci.state.max_abs_K : 1e100;
    double rmsK = isfinite(out.ricci.state.rms_K) ? out.ricci.state.rms_K : 1e100;
    double spread = isfinite(out.develop.global_vertex_max_spread) ? out.develop.global_vertex_max_spread : 1e100;
    double period = isfinite(out.develop.period_fit_rms) ? out.develop.period_fit_rms : 1e100;
    double overlap = isfinite(out.packing.nonedge_max_overlap) ? out.packing.nonedge_max_overlap : 1e100;
    double frac = isfinite(out.packing.nonedge_overlap_fraction) ? out.packing.nonedge_overlap_fraction : 1.0;
    double ratio = isfinite(out.ricci.radius_ratio) && out.ricci.radius_ratio > 0.0 ? out.ricci.radius_ratio : 1e100;
    double delta = isfinite(out.packing.radius_stats.polydispersity_delta) ? out.packing.radius_stats.polydispersity_delta : 1e100;
    auto logterm = [](double x, double scale) {
        scale = max(scale, 1e-300);
        return log10(1.0 + max(0.0, x) / scale);
    };

    double radius_penalty = pow(max(0.0, log(max(ratio / cfg.target_max_ratio, 1e-300))), 2.0);
    double delta_penalty = pow(max(0.0, delta - cfg.target_delta), 2.0);

    ScoreBreakdown s;
    s.K = cfg.wK * logterm(maxK, cfg.K_scale);
    if (cfg.surgery_metric == "v6") {
        s.rmsK = cfg.wRmsK * logterm(rmsK, cfg.rmsK_scale);
        s.develop = cfg.wDev * logterm(spread, cfg.dev_scale);
        s.period = cfg.wPer * logterm(period, cfg.period_scale);
        s.overlap = cfg.wOv * logterm(overlap, cfg.overlap_scale);
        s.overlap_count = cfg.wOvN * frac;
        s.crossing = cfg.wCross * static_cast<double>(out.packing.edge_crossing_count);
        s.degenerate = cfg.wDeg * static_cast<double>(out.packing.face_degenerate_count);
        s.orientation = cfg.wOrient * static_cast<double>(out.packing.orientation_neighbor_violation_count);
        s.sort = cfg.wSort * (isfinite(out.dist.D_sort) ? out.dist.D_sort : 1e100);
        s.delta = cfg.wDelta * fabs(delta - cfg.target_delta);
        s.ratio = cfg.wRatio * max(0.0, log(max(ratio / cfg.target_radius_ratio, 1e-300)));
    } else {
        s.develop = cfg.wDev * logterm(spread, cfg.dev_scale);
        s.period = cfg.wPer * logterm(period, cfg.period_scale);
        s.overlap = cfg.wOv * logterm(overlap, cfg.overlap_scale);
        s.overlap_count = cfg.wOvN * frac;
        s.crossing = cfg.wCross * static_cast<double>(out.packing.edge_crossing_count);
        s.degenerate = cfg.wDeg * static_cast<double>(out.packing.face_degenerate_count);
        s.orientation = cfg.wOrient * static_cast<double>(out.packing.orientation_neighbor_violation_count);
        s.radius = cfg.wRad * radius_penalty;
        s.delta = cfg.wDelta * delta_penalty;
    }
    s.total = s.K + s.develop + s.period + s.overlap + s.overlap_count +
              s.crossing + s.degenerate + s.orientation + s.radius + s.delta +
              s.rmsK + s.sort + s.ratio;
    return s;
}

static double surgery_score(const RunOutput& out, const Config& cfg) {
    return compute_score_breakdown(out, cfg).total;
}

static RunOutput relax_mesh_from_u(Config cfg, Mesh mesh, vector<double> u, const vector<double>& u0,
                                   bool keep_stats, bool verbose) {
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
        out.ricci = run_ricci(out.mesh, std::move(u), cfg, u0, keep_stats, verbose);
        out.time.ricci = now_sec() - t0;
        out.dist = compute_distribution_diagnostics(out.ricci.state.u, u0, cfg,
                                                    out.ricci.state.max_abs_K,
                                                    out.ricci.state.rms_K);
        if (!u0.empty()) {
            vector<int> target = assign_target_degrees(u0, cfg);
            out.degree = compute_degree_diagnostics(out.mesh, u0, target, cfg);
        }
        if ((cfg.surgery_metric == "linear" || cfg.linear_pareto_scan || cfg.screen_topologies) && !out.dist.u0.empty()) {
            out.linear = compute_linear_response(out.mesh, out.dist.u0, cfg);
        }
        if (cfg.develop) {
            t0 = now_sec();
            out.develop = develop_mesh(out.mesh, out.ricci.state.r);
            out.time.develop = now_sec() - t0;
            out.packing = compute_packing_diagnostics(out.mesh, out.ricci, out.develop, cfg, out.dist);
        }
        out.time.total = now_sec() - t_total0;
    } catch (const exception& e) {
        out.valid = false;
        out.error = e.what();
        out.time.total = now_sec() - t_total0;
    }
    return out;
}

static RunOutput evaluate_linear_mesh(Config cfg, Mesh mesh, const vector<double>& u0) {
    cfg.radius_mode = "fixed";
    cfg.surgery_metric = "linear";
    cfg.max_iter = 0;
    cfg.develop = false;
    return relax_mesh_from_u(cfg, std::move(mesh), u0, u0, false, false);
}

static RunOutput refine_after_linear_search(const RunOutput& lin_out, Config cfg) {
    RunOutput out;
    out.cfg = cfg;
    out.mesh = lin_out.mesh;
    out.flips_accepted = lin_out.flips_accepted;
    out.linear = lin_out.linear;
    out.degree = lin_out.degree;
    double t_total0 = now_sec();
    try {
        if (lin_out.dist.u0.empty() || lin_out.linear.du.empty()) throw runtime_error("missing linear predictor for refinement");
        vector<double> u = lin_out.dist.u0;
        double max_initial_du = 0.0;
        for (double x : lin_out.linear.du) max_initial_du = max(max_initial_du, fabs(x));
        double scale = max_initial_du > cfg.refine_du_bound
                       ? cfg.refine_du_bound / max(max_initial_du, 1e-300)
                       : 1.0;
        for (int i = 0; i < static_cast<int>(u.size()); ++i) u[i] += scale * lin_out.linear.du[i];
        normalize_u(u);

        State current, candidate;
        compute_state_inplace(out.mesh, u, current);
        double dt = 0.01;
        vector<double> cand(u.size(), 0.0);
        int boundary_rejections = 0;
        int iter = 0;
        for (iter = 1; iter <= cfg.refine_max_iter; ++iter) {
            bool accepted = false;
            while (!accepted) {
                if (dt < 1e-14) {
                    iter--;
                    goto done_refine;
                }
                for (int i = 0; i < static_cast<int>(u.size()); ++i) cand[i] = u[i] - dt * current.K[i];
                normalize_u(cand);
                bool inside = true;
                for (int i = 0; i < static_cast<int>(u.size()); ++i) {
                    if (fabs(cand[i] - lin_out.dist.u0[i]) > cfg.refine_du_bound) {
                        inside = false;
                        break;
                    }
                }
                if (!inside) {
                    boundary_rejections++;
                    dt *= 0.5;
                    continue;
                }
                compute_state_inplace(out.mesh, cand, candidate);
                if (candidate.E_K < current.E_K) {
                    u.swap(cand);
                    current = std::move(candidate);
                    dt = min(dt * 1.05, 1.0);
                    accepted = true;
                } else {
                    dt *= 0.5;
                }
            }
            if (current.max_abs_K < cfg.refine_tol) break;
        }
done_refine:
        out.ricci.state = std::move(current);
        out.ricci.iterations = max(0, iter);
        out.ricci.boundary_rejections = boundary_rejections;
        out.ricci.stopped_reason = out.ricci.state.max_abs_K < cfg.refine_tol
                                    ? "linear refinement converged"
                                    : "linear refinement stopped";
        finalize_ricci(out.ricci);
        out.ricci.final_E_reg = out.ricci.state.E_K;
        out.dist = compute_distribution_diagnostics(out.ricci.state.u, lin_out.dist.u0, cfg,
                                                    out.ricci.state.max_abs_K,
                                                    out.ricci.state.rms_K);
        if (!lin_out.degree.target_degree.empty()) {
            out.degree = compute_degree_diagnostics(out.mesh, lin_out.dist.u0, lin_out.degree.target_degree, cfg);
        }
        if (cfg.develop) {
            out.develop = develop_mesh(out.mesh, out.ricci.state.r);
            out.packing = compute_packing_diagnostics(out.mesh, out.ricci, out.develop, cfg, out.dist);
            double max_du = 0.0;
            for (int i = 0; i < static_cast<int>(out.ricci.state.u.size()); ++i) {
                max_du = max(max_du, fabs(out.ricci.state.u[i] - lin_out.dist.u0[i]));
            }
            out.packing.physical_INC_like_candidate_v6 = out.packing.physical_INC_like_candidate_v6 &&
                                                         out.dist.D_sort < cfg.dist_tol_sort &&
                                                         max_du <= cfg.refine_du_bound + 1e-12;
        }
        out.time.total = now_sec() - t_total0;
    } catch (const exception& e) {
        out.valid = false;
        out.error = e.what();
        out.time.total = now_sec() - t_total0;
    }
    return out;
}

static RunOutput run_linear_surgery(Config cfg) {
    double t0 = now_sec();
    filesystem::create_directories(cfg.out);
    Config base = cfg;
    base.radius_mode = "fixed";
    base.surgery_metric = "linear";
    base.max_iter = 0;
    base.develop = false;
    RunOutput current = execute_run(base, false, false);
    if (!current.valid) return current;
    int initial_flips_accepted = current.flips_accepted;
    int surgery_accepts = 0;
    double current_score = linear_total_score(current, cfg);

    ofstream stats(filesystem::path(cfg.out) / "linear_surgery_stats.csv");
    stats << "step,accepted,S_old,S_new,rms_K0,max_K0,rms_K_lin,max_K_lin,du_rms,du_max,"
             "degree_radius_corr,singular_response_flag,flipped_edge_i,flipped_edge_j,time_step_sec\n";
    stats << scientific << setprecision(17);

    mt19937_64 rng(cfg.seed + 0x51ed1234ULL);
    uniform_real_distribution<double> uni(0.0, 1.0);
    for (int step = 1; step <= cfg.linear_surgery_steps; ++step) {
        double ts = now_sec();
        double old_score = current_score;
        RunOutput best;
        double best_score = numeric_limits<double>::infinity();
        bool have_best = false;
        int best_i = -1, best_j = -1;

        for (int trial = 0; trial < cfg.linear_surgery_trials_per_step; ++trial) {
            Mesh proposal_mesh = current.mesh;
            int ei = -1, ej = -1;
            bool flipped = attempt_linear_guided_flip(proposal_mesh, rng, cfg, current, ei, ej);
            if (!flipped) continue;
            string error;
            if (!validate_mesh(proposal_mesh, &error)) continue;
            RunOutput proposal = evaluate_linear_mesh(cfg, std::move(proposal_mesh), current.dist.u0);
            if (!proposal.valid || !proposal.linear.ran) continue;
            double S = linear_total_score(proposal, cfg);
            if (S < best_score) {
                best_score = S;
                best = std::move(proposal);
                best_i = ei;
                best_j = ej;
                have_best = true;
            }
        }

        bool accepted = false;
        if (have_best && best_score < current_score) {
            accepted = true;
        } else if (have_best && cfg.linear_surgery_temperature > 0.0) {
            double p = exp(-(best_score - current_score) / cfg.linear_surgery_temperature);
            accepted = uni(rng) < p;
        }
        int flipped_i = -1, flipped_j = -1;
        if (accepted) {
            current = std::move(best);
            current_score = best_score;
            flipped_i = best_i;
            flipped_j = best_j;
            surgery_accepts++;
        }

        stats << step << ',' << (accepted ? 1 : 0) << ',' << old_score << ',' << current_score << ','
              << current.linear.rms_K0 << ',' << current.linear.max_K0 << ','
              << current.linear.rms_K_lin << ',' << current.linear.max_K_lin << ','
              << current.linear.du_rms << ',' << current.linear.du_max << ','
              << current.linear.degree_radius_corr << ','
              << (current.linear.singular_response_flag ? 1 : 0) << ','
              << flipped_i << ',' << flipped_j << ',' << (now_sec() - ts) << '\n';
    }
    current.cfg = cfg;
    current.flips_accepted = initial_flips_accepted + surgery_accepts;
    current.time.total = now_sec() - t0;
    if (cfg.refine_after_linear) {
        RunOutput refined = refine_after_linear_search(current, cfg);
        refined.cfg = cfg;
        refined.flips_accepted = current.flips_accepted;
        refined.linear = current.linear;
        refined.time.total = now_sec() - t0;
        return refined;
    }
    return current;
}

static RunOutput run_surgery(Config cfg) {
    double t_surgery0 = now_sec();
    filesystem::create_directories(cfg.out);
    cfg.develop = true;
    bool quiet_requested = cfg.progress_every == 0;
    cfg.progress_every = 0;
    if (quiet_requested) cfg.surgery_report_every = 0;
    Config base_cfg = cfg;
    RunOutput current = execute_run(base_cfg, false, false);
    if (!current.valid) {
        current.time.total = now_sec() - t_surgery0;
        return current;
    }
    int initial_flips_accepted = current.flips_accepted;
    int surgery_accepts = 0;
    double current_score = surgery_score(current, cfg);

    ofstream stats(filesystem::path(cfg.out) / "surgery_stats.csv");
    stats << "step,accepted,score_old,score_new,max_abs_K,develop_spread,period_fit_rms,"
             "nonedge_max_overlap,nonedge_overlap_count,edge_crossing_count,face_degenerate_count,"
             "face_negative_signed_count,orientation_neighbor_violation_count,radius_ratio,"
             "polydispersity_delta,degree6_fraction,physical_INC_like_candidate,time_step_sec,"
             "radius_mode,D_u,D_sort,D_quantile,D_delta,D_ratio,target_delta,target_radius_ratio,"
             "fixed_intrinsic_good,weak_distribution_preserved,physical_INC_like_candidate_v6\n";
    stats << scientific << setprecision(17);

    mt19937_64 rng(cfg.seed + 0x9e3779b97f4a7c15ULL);
    uniform_real_distribution<double> uni(0.0, 1.0);
    for (int step = 1; step <= cfg.surgery_steps; ++step) {
        double t_step = now_sec();
        double old_score = current_score;
        RunOutput best;
        double best_score = numeric_limits<double>::infinity();
        bool have_best = false;

        for (int trial = 0; trial < cfg.surgery_trials_per_step; ++trial) {
            Mesh proposal_mesh = current.mesh;
            bool flipped = false;
            for (int bf = 0; bf < cfg.surgery_batch_flips; ++bf) {
                bool ok = (bf == 0) ? attempt_biased_flip(proposal_mesh, rng, cfg, current)
                                    : attempt_random_flip(proposal_mesh, rng);
                flipped = flipped || ok;
            }
            if (!flipped) continue;
            string error;
            if (!validate_mesh(proposal_mesh, &error)) continue;

            Config trial_cfg = cfg;
            trial_cfg.max_iter = min(cfg.max_iter, cfg.surgery_ricci_iter);
            trial_cfg.progress_every = 0;
            RunOutput proposal = relax_mesh_from_u(trial_cfg, std::move(proposal_mesh),
                                                   current.ricci.state.u, current.dist.u0,
                                                   false, false);
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

        double step_sec = now_sec() - t_step;
        stats << step << ',' << (accepted ? 1 : 0) << ',' << old_score << ',' << current_score << ','
              << current.ricci.state.max_abs_K << ',' << current.develop.global_vertex_max_spread << ','
              << current.develop.period_fit_rms << ',' << current.packing.nonedge_max_overlap << ','
              << current.packing.nonedge_overlap_count << ',' << current.packing.edge_crossing_count << ','
              << current.packing.face_degenerate_count << ',' << current.packing.face_negative_signed_count << ','
              << current.packing.orientation_neighbor_violation_count << ',' << current.ricci.radius_ratio << ','
              << current.packing.radius_stats.polydispersity_delta << ',' << current.packing.degree6_fraction << ','
              << (current.packing.physical_INC_like_candidate ? 1 : 0) << ',' << step_sec << ','
              << cfg.radius_mode << ',' << current.dist.D_u << ',' << current.dist.D_sort << ','
              << current.dist.D_quantile << ',' << current.dist.D_delta << ',' << current.dist.D_ratio << ','
              << cfg.target_delta << ',' << cfg.target_radius_ratio << ','
              << (current.packing.fixed_intrinsic_good ? 1 : 0) << ','
              << (current.packing.weak_distribution_preserved ? 1 : 0) << ','
              << (current.packing.physical_INC_like_candidate_v6 ? 1 : 0) << '\n';

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
    {
        ScoreBreakdown sb = compute_score_breakdown(current, cfg);
        ofstream f(filesystem::path(cfg.out) / "best_candidate_summary.txt");
        f << scientific << setprecision(17);
        f << "final_score: " << sb.total << "\n";
        f << "score_K: " << sb.K << "\n";
        f << "score_rmsK: " << sb.rmsK << "\n";
        f << "score_develop: " << sb.develop << "\n";
        f << "score_period: " << sb.period << "\n";
        f << "score_overlap: " << sb.overlap << "\n";
        f << "score_overlap_count: " << sb.overlap_count << "\n";
        f << "score_crossing: " << sb.crossing << "\n";
        f << "score_degenerate: " << sb.degenerate << "\n";
        f << "score_orientation: " << sb.orientation << "\n";
        f << "score_radius: " << sb.radius << "\n";
        f << "score_delta: " << sb.delta << "\n";
        f << "score_sort: " << sb.sort << "\n";
        f << "score_ratio: " << sb.ratio << "\n";
        f << "intrinsic_endpoint: " << (current.packing.intrinsic_endpoint ? "yes" : "no") << "\n";
        f << "developed_consistent: " << (current.packing.developed_consistent ? "yes" : "no") << "\n";
        f << "local_theta_valid: " << (current.packing.local_theta_valid ? "yes" : "no") << "\n";
        f << "contact_edges_valid: " << (current.packing.contact_edges_valid ? "yes" : "no") << "\n";
        f << "nonedge_overlap_free: " << (current.packing.nonedge_overlap_free ? "yes" : "no") << "\n";
        f << "geometry_valid: " << (current.packing.topology_geometry_valid ? "yes" : "no") << "\n";
        f << "physically_reasonable_radius: " << (current.packing.physically_reasonable_radius ? "yes" : "no") << "\n";
        f << "valid_disk_packing_candidate: " << (current.packing.valid_disk_packing_candidate ? "yes" : "no") << "\n";
        f << "strong_INC_candidate: " << (current.packing.strong_INC_candidate ? "yes" : "no") << "\n";
        f << "physical_INC_like_candidate: " << (current.packing.physical_INC_like_candidate ? "yes" : "no") << "\n";
        f << "fixed_intrinsic_good: " << (current.packing.fixed_intrinsic_good ? "yes" : "no") << "\n";
        f << "weak_distribution_preserved: " << (current.packing.weak_distribution_preserved ? "yes" : "no") << "\n";
        f << "distribution_preserved: " << (current.packing.distribution_preserved ? "yes" : "no") << "\n";
        f << "physical_INC_like_candidate_v6: " << (current.packing.physical_INC_like_candidate_v6 ? "yes" : "no") << "\n";
        f << "radius_ratio: " << current.ricci.radius_ratio << "\n";
        f << "polydispersity_delta: " << current.packing.radius_stats.polydispersity_delta << "\n";
        f << "\nDistribution preservation diagnostics\n";
        f << "radius_mode: " << cfg.radius_mode << "\n";
        f << "target_delta: " << cfg.target_delta << "\n";
        f << "target_radius_ratio: " << cfg.target_radius_ratio << "\n";
        f << "D_u: " << current.dist.D_u << "\n";
        f << "D_sort: " << current.dist.D_sort << "\n";
        f << "D_quantile: " << current.dist.D_quantile << "\n";
        f << "D_delta: " << current.dist.D_delta << "\n";
        f << "D_ratio: " << current.dist.D_ratio << "\n";
        f << "nonedge_overlap_count: " << current.packing.nonedge_overlap_count << "\n";
        f << "nonedge_max_overlap: " << current.packing.nonedge_max_overlap << "\n";
        f << "edge_crossing_count: " << current.packing.edge_crossing_count << "\n";
        f << "face_degenerate_count: " << current.packing.face_degenerate_count << "\n";
        f << "face_negative_signed_count: " << current.packing.face_negative_signed_count << "\n";
        f << "orientation_neighbor_violation_count: " << current.packing.orientation_neighbor_violation_count << "\n";
        f << "degree6_fraction: " << current.packing.degree6_fraction << "\n";
        f << "noncrystalline_topology_heuristic: " << (current.packing.noncrystalline_topology_heuristic ? "yes" : "no") << "\n";
        f << "\nInterpretation relative to NM/PRL\n";
        if (!current.packing.physical_INC_like_candidate_v6) {
            f << "recommended_interpretation: mathematical circle-packing candidate, not physical INC-like unless all v6 flags are satisfied.\n";
        } else {
            f << "recommended_interpretation: v6 physical INC-like candidate under the implemented combinatorial/development diagnostics.\n";
        }
    }
    current.time.total = now_sec() - t_surgery0;
    return current;
}

static vector<RunOutput> run_multistart_surgery(Config cfg) {
    vector<RunOutput> results;
    results.reserve(cfg.multistart_count);
    ofstream f(cfg.multistart_out);
    f << "run_id,seed,flips_requested,flips_accepted,final_score,final_max_abs_K,develop_spread,"
         "period_fit_rms,nonedge_overlap_count,nonedge_max_overlap,edge_crossing_count,"
         "face_degenerate_count,orientation_neighbor_violation_count,radius_ratio,polydispersity_delta,"
         "valid_disk_packing_candidate,strong_INC_candidate,physical_INC_like_candidate,"
         "physical_INC_like_candidate_v6,time_total_sec\n";
    f << scientific << setprecision(17);

    for (int run = 0; run < cfg.multistart_count; ++run) {
        Config local = cfg;
        local.seed = cfg.seed + static_cast<unsigned long long>(run);
        local.surgery = true;
        local.multistart_surgery = false;
        local.out = cfg.out + "_multistart_" + to_string(run);
        RunOutput out = run_surgery(local);
        ScoreBreakdown sb = compute_score_breakdown(out, local);
        f << run << ',' << local.seed << ',' << local.flips << ',' << out.flips_accepted << ','
          << sb.total << ',' << out.ricci.state.max_abs_K << ',' << out.develop.global_vertex_max_spread << ','
          << out.develop.period_fit_rms << ',' << out.packing.nonedge_overlap_count << ','
          << out.packing.nonedge_max_overlap << ',' << out.packing.edge_crossing_count << ','
          << out.packing.face_degenerate_count << ',' << out.packing.orientation_neighbor_violation_count << ','
          << out.ricci.radius_ratio << ',' << out.packing.radius_stats.polydispersity_delta << ','
          << (out.packing.valid_disk_packing_candidate ? 1 : 0) << ','
          << (out.packing.strong_INC_candidate ? 1 : 0) << ','
          << (out.packing.physical_INC_like_candidate ? 1 : 0) << ','
          << (out.packing.physical_INC_like_candidate_v6 ? 1 : 0) << ','
          << out.time.total << '\n';
        if (cfg.progress_every > 0) {
            cout << scientific << setprecision(6)
                 << "multistart run=" << run
                 << " seed=" << local.seed
                 << " score=" << sb.total
                 << " maxK=" << out.ricci.state.max_abs_K
                 << " ratio=" << out.ricci.radius_ratio
                 << " physical=" << (out.packing.physical_INC_like_candidate ? 1 : 0) << '\n';
        }
        results.push_back(std::move(out));
    }
    return results;
}

static vector<RunOutput> run_pareto_scan(Config cfg) {
    vector<RunOutput> rows;
    ofstream f(cfg.pareto_out);
    f << "mode,target_delta,target_radius_ratio,seed,flips_requested,flips_accepted,final_score,"
         "final_max_abs_K,final_rms_K,develop_spread,period_fit_rms,nonedge_overlap_count,"
         "nonedge_max_overlap,edge_crossing_count,face_degenerate_count,orientation_neighbor_violation_count,"
         "radius_ratio,polydispersity_delta,D_sort,D_quantile,valid_disk_packing_candidate,"
         "strong_INC_candidate,physical_INC_like_candidate,fixed_intrinsic_good,weak_distribution_preserved,"
         "physical_INC_like_candidate_v6,time_total_sec\n";
    f << scientific << setprecision(17);

    for (const string& mode : cfg.pareto_modes) {
        for (double delta : cfg.pareto_delta_list) {
            for (double ratio : cfg.pareto_ratio_list) {
                for (unsigned long long seed : cfg.pareto_seeds) {
                    Config local = cfg;
                    local.radius_mode = mode;
                    local.target_delta = delta;
                    local.target_radius_ratio = ratio;
                    local.seed = seed;
                    local.radius_seed = cfg.radius_seed + 1000003ULL * seed +
                                        static_cast<unsigned long long>(llround(delta * 1000.0)) * 9176ULL +
                                        static_cast<unsigned long long>(llround(ratio * 100.0));
                    local.progress_every = 0;
                    local.scan_write_details = cfg.scan_write_details;
                    local.out = cfg.out + "_pareto_" + mode + "_d" + to_string(static_cast<int>(round(delta * 100))) +
                                "_r" + to_string(static_cast<int>(round(ratio))) +
                                "_s" + to_string(seed);
                    apply_output_defaults(local, true);
                    RunOutput out = local.surgery ? run_surgery(local) : execute_run(local, cfg.scan_write_details, false);
                    ScoreBreakdown sb = compute_score_breakdown(out, local);
                    f << mode << ',' << delta << ',' << ratio << ',' << seed << ','
                      << local.flips << ',' << out.flips_accepted << ',' << sb.total << ','
                      << out.ricci.state.max_abs_K << ',' << out.ricci.state.rms_K << ','
                      << out.develop.global_vertex_max_spread << ',' << out.develop.period_fit_rms << ','
                      << out.packing.nonedge_overlap_count << ',' << out.packing.nonedge_max_overlap << ','
                      << out.packing.edge_crossing_count << ',' << out.packing.face_degenerate_count << ','
                      << out.packing.orientation_neighbor_violation_count << ',' << out.ricci.radius_ratio << ','
                      << out.packing.radius_stats.polydispersity_delta << ',' << out.dist.D_sort << ','
                      << out.dist.D_quantile << ','
                      << (out.packing.valid_disk_packing_candidate ? 1 : 0) << ','
                      << (out.packing.strong_INC_candidate ? 1 : 0) << ','
                      << (out.packing.physical_INC_like_candidate ? 1 : 0) << ','
                      << (out.packing.fixed_intrinsic_good ? 1 : 0) << ','
                      << (out.packing.weak_distribution_preserved ? 1 : 0) << ','
                      << (out.packing.physical_INC_like_candidate_v6 ? 1 : 0) << ','
                      << out.time.total << '\n';
                    if (cfg.scan_write_details) write_outputs(out);
                    rows.push_back(std::move(out));
                }
            }
        }
    }
    return rows;
}

static void run_topology_screen(Config cfg) {
    ofstream f(cfg.screen_out);
    f << "sample,seed,flips_requested,flips_accepted,S_pred,rms_K0,max_K0,rms_K_lin,max_K_lin,"
         "du_rms,du_max,response_ratio,degree_radius_corr,singular_response_flag,noncrystalline_topology_heuristic\n";
    f << scientific << setprecision(17);

    vector<double> shared_u0;
    if (cfg.screen_same_radii) shared_u0 = generate_target_log_radii(cfg.nx * cfg.ny, cfg);
    for (int sample = 0; sample < cfg.screen_count; ++sample) {
        unsigned long long seed = cfg.seed + static_cast<unsigned long long>(sample);
        try {
            mt19937_64 rng(seed);
            Mesh mesh = generate_initial_mesh(cfg.nx, cfg.ny);
            int accepted = randomize_topology(mesh, cfg.screen_flips, cfg.validate_every, rng);
            string error;
            if (!validate_mesh(mesh, &error)) throw runtime_error(error);
            vector<double> u0 = cfg.screen_same_radii ? shared_u0 : generate_target_log_radii(mesh.N, cfg);
            LinearResponseResult lr = compute_linear_response(mesh, u0, cfg);
            vector<int> degree = vertex_degrees(mesh);
            int degree6 = 0;
            for (int d : degree) if (d == 6) degree6++;
            bool noncrystalline = (static_cast<double>(degree6) / mesh.N) < 0.95;
            f << sample << ',' << seed << ',' << cfg.screen_flips << ',' << accepted << ','
              << lr.S_pred << ',' << lr.rms_K0 << ',' << lr.max_K0 << ','
              << lr.rms_K_lin << ',' << lr.max_K_lin << ',' << lr.du_rms << ','
              << lr.du_max << ',' << lr.response_ratio << ',' << lr.degree_radius_corr << ','
              << (lr.singular_response_flag ? 1 : 0) << ',' << (noncrystalline ? 1 : 0) << '\n';
        } catch (...) {
            f << sample << ',' << seed << ',' << cfg.screen_flips << ",0,"
              << "1e100,1e100,1e100,1e100,1e100,1e100,1e100,1e100,0,1,0\n";
        }
    }
}

static vector<RunOutput> run_linear_pareto_scan(Config cfg) {
    vector<RunOutput> rows;
    ofstream f(cfg.pareto_out);
    f << "target_delta,target_radius_ratio,seed,S_pred,rms_K_lin,max_K_lin,du_rms,du_max,"
         "degree_radius_corr,linear_predicted_good,refined_physical_candidate\n";
    f << scientific << setprecision(17);

    for (double delta : cfg.pareto_delta_list) {
        for (double ratio : cfg.pareto_ratio_list) {
            for (unsigned long long seed : cfg.pareto_seeds) {
                Config local = cfg;
                local.target_delta = delta;
                local.target_radius_ratio = ratio;
                local.seed = seed;
                local.radius_seed = cfg.radius_seed + 1000003ULL * seed +
                                    static_cast<unsigned long long>(llround(delta * 1000.0)) * 9176ULL +
                                    static_cast<unsigned long long>(llround(ratio * 100.0));
                local.radius_mode = "fixed";
                local.surgery_metric = "linear";
                local.out = cfg.out + "_linear_pareto_d" + to_string(static_cast<int>(round(delta * 100))) +
                            "_r" + to_string(static_cast<int>(round(ratio))) +
                            "_s" + to_string(seed);
                RunOutput out;
                if (local.linear_surgery_steps > 0) {
                    out = run_linear_surgery(local);
                } else {
                    out = execute_run(local, false, false);
                }
                bool refined_physical = out.packing.physical_INC_like_candidate_v6;
                double S = out.linear.ran ? out.linear.S_pred : 1e100;
                f << delta << ',' << ratio << ',' << seed << ',' << S << ','
                  << (out.linear.ran ? out.linear.rms_K_lin : 1e100) << ','
                  << (out.linear.ran ? out.linear.max_K_lin : 1e100) << ','
                  << (out.linear.ran ? out.linear.du_rms : 1e100) << ','
                  << (out.linear.ran ? out.linear.du_max : 1e100) << ','
                  << (out.linear.ran ? out.linear.degree_radius_corr : 0.0) << ','
                  << (out.linear.linear_predicted_good ? 1 : 0) << ','
                  << (refined_physical ? 1 : 0) << '\n';
                rows.push_back(std::move(out));
            }
        }
    }
    return rows;
}

struct DegreeProfile {
    string name;
    double f4, f5, f6, f7, f8;
};

static vector<DegreeProfile> degree_profiles() {
    return {
        {"weak", 0.04, 0.16, 0.60, 0.16, 0.04},
        {"medium", 0.08, 0.17, 0.50, 0.17, 0.08},
        {"strong", 0.12, 0.18, 0.40, 0.18, 0.12},
        {"very_strong", 0.16, 0.18, 0.32, 0.18, 0.16}
    };
}

static RunOutput run_v8_design_pipeline(Config cfg) {
    cfg.radius_mode = "fixed";
    cfg.surgery_metric = "linear";
    if (cfg.degree_target_mode == "none") cfg.degree_target_mode = "quantile";
    if (cfg.proposal_mode == "random") cfg.proposal_mode = "mismatch";
    if (cfg.linear_surgery_steps > 0) return run_linear_surgery(cfg);
    return execute_run(cfg, false, false);
}

static vector<RunOutput> run_degree_fraction_scan(Config cfg) {
    vector<RunOutput> rows;
    ofstream f(cfg.degree_fraction_out);
    f << "profile,f4,f5,f6,f7,f8,S_degree,S_pred,S_total,rms_K_lin,max_K_lin,du_rms,du_max,"
         "degree_match_fraction,degree_radius_corr,linear_predicted_good,refined_physical_candidate\n";
    f << scientific << setprecision(17);

    for (const DegreeProfile& p : degree_profiles()) {
        Config local = cfg;
        local.radius_mode = "fixed";
        local.surgery_metric = "linear";
        local.degree_target_mode = "quantile";
        local.frac_deg4 = p.f4; local.frac_deg5 = p.f5; local.frac_deg6 = p.f6;
        local.frac_deg7 = p.f7; local.frac_deg8 = p.f8;
        if (local.proposal_mode == "random") local.proposal_mode = "mismatch";
        local.out = cfg.out + "_degree_profile_" + p.name;
        RunOutput out = run_v8_design_pipeline(local);
        double total = linear_total_score(out, local);
        f << p.name << ',' << p.f4 << ',' << p.f5 << ',' << p.f6 << ',' << p.f7 << ',' << p.f8 << ','
          << out.degree.S_degree << ',' << out.linear.S_pred << ',' << total << ','
          << out.linear.rms_K_lin << ',' << out.linear.max_K_lin << ','
          << out.linear.du_rms << ',' << out.linear.du_max << ','
          << out.degree.degree_match_fraction << ',' << out.degree.degree_radius_corr << ','
          << (out.linear.linear_predicted_good ? 1 : 0) << ','
          << (out.packing.physical_INC_like_candidate_v6 ? 1 : 0) << '\n';
        rows.push_back(std::move(out));
    }
    return rows;
}

static vector<RunOutput> run_strategy_compare(Config cfg) {
    vector<pair<string, Config>> tasks;
    Config base = cfg;
    base.radius_mode = "fixed";
    base.surgery_metric = "linear";
    base.degree_target_mode = "quantile";
    if (base.proposal_mode == "random") base.proposal_mode = "mismatch";

    Config random_only = base;
    random_only.degree_precondition = false;
    random_only.linear_surgery_steps = 0;
    tasks.push_back({"random_flips_only", random_only});

    Config degree_only = base;
    degree_only.degree_precondition = true;
    degree_only.linear_surgery_steps = 0;
    tasks.push_back({"degree_precondition_only", degree_only});

    Config degree_linear = base;
    degree_linear.degree_precondition = true;
    tasks.push_back({"degree_precondition_plus_linear_surgery", degree_linear});

    Config linear_only = base;
    linear_only.degree_precondition = false;
    tasks.push_back({"linear_surgery_only", linear_only});

    vector<RunOutput> rows;
    ofstream f(cfg.strategy_out);
    f << "strategy,S_degree,S_pred,S_total,du_rms,du_max,rms_K_lin,max_K_lin,degree_radius_corr,"
         "degree_match_fraction,linear_predicted_good,refined_physical_candidate\n";
    f << scientific << setprecision(17);
    for (auto& item : tasks) {
        Config local = item.second;
        local.out = cfg.out + "_strategy_" + item.first;
        RunOutput out = run_v8_design_pipeline(local);
        double total = linear_total_score(out, local);
        f << item.first << ',' << out.degree.S_degree << ',' << out.linear.S_pred << ',' << total << ','
          << out.linear.du_rms << ',' << out.linear.du_max << ','
          << out.linear.rms_K_lin << ',' << out.linear.max_K_lin << ','
          << out.degree.degree_radius_corr << ',' << out.degree.degree_match_fraction << ','
          << (out.linear.linear_predicted_good ? 1 : 0) << ','
          << (out.packing.physical_INC_like_candidate_v6 ? 1 : 0) << '\n';
        rows.push_back(std::move(out));
    }
    return rows;
}

static bool run_tests() {
    bool ok = true;
    cout << "Running v8 tests\n";
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 0; cfg.seed = 1; cfg.radius_seed = 1;
        cfg.radius_mode = "fixed"; cfg.target_delta = 0.0; cfg.target_radius_ratio = 5.0;
        cfg.surgery_metric = "linear"; cfg.degree_target_mode = "none"; cfg.progress_every = 0;
        RunOutput out = execute_run(cfg, false, false);
        bool pass = out.valid && out.linear.ran &&
                    out.linear.rms_K0 < 1e-12 &&
                    out.linear.du_max < 1e-10 &&
                    out.linear.S_pred < 1e-8 &&
                    out.linear.linear_predicted_good;
        ok = ok && pass;
        cout << "Test 1 regular equal no degree target: " << (pass ? "PASS" : "FAIL")
             << " rmsK0=" << scientific << setprecision(6) << out.linear.rms_K0
             << " du_max=" << out.linear.du_max
             << " S=" << out.linear.S_pred << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 0; cfg.seed = 2; cfg.radius_seed = 101;
        cfg.radius_mode = "fixed"; cfg.target_dist = "lognormal"; cfg.target_delta = 0.25;
        cfg.target_radius_ratio = 5.0; cfg.surgery_metric = "linear";
        cfg.degree_target_mode = "quantile"; cfg.progress_every = 0;
        RunOutput out = execute_run(cfg, false, false);
        int sum_target = accumulate(out.degree.target_degree.begin(), out.degree.target_degree.end(), 0);
        bool pass = out.valid && out.degree.ran &&
                    sum_target == 6 * out.mesh.N &&
                    out.degree.target_degree_radius_corr > 0.5;
        ok = ok && pass;
        cout << "Test 2 target degree assignment: " << (pass ? "PASS" : "FAIL")
             << " sum_target=" << sum_target
             << " target_corr=" << out.degree.target_degree_radius_corr << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 100; cfg.seed = 1; cfg.radius_seed = 101;
        cfg.radius_mode = "fixed"; cfg.target_dist = "lognormal"; cfg.target_delta = 0.25;
        cfg.target_radius_ratio = 5.0; cfg.surgery_metric = "linear"; cfg.degree_target_mode = "quantile";
        cfg.progress_every = 0;
        RunOutput initial = execute_run(cfg, false, false);
        cfg.degree_precondition = true; cfg.degree_pre_steps = 50; cfg.degree_pre_trials = 20;
        cfg.out = "degree_pre_v8_test";
        RunOutput final = execute_run(cfg, false, false);
        bool pass = initial.valid && final.valid && final.degree.S_degree <= initial.degree.S_degree + 1e-12;
        ok = ok && pass;
        cout << "Test 3 degree precondition small: " << (pass ? "PASS" : "FAIL")
             << " S_initial=" << initial.degree.S_degree
             << " S_final=" << final.degree.S_degree << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 100; cfg.seed = 1; cfg.radius_seed = 101;
        cfg.radius_mode = "fixed"; cfg.target_dist = "lognormal"; cfg.target_delta = 0.25;
        cfg.target_radius_ratio = 5.0; cfg.surgery_metric = "linear"; cfg.degree_target_mode = "quantile";
        cfg.degree_precondition = true; cfg.degree_pre_steps = 50; cfg.degree_pre_trials = 20;
        cfg.out = "degree_linear_v8_test"; cfg.progress_every = 0;
        RunOutput out = execute_run(cfg, true, false);
        write_outputs(out);
        bool pass = out.valid && out.linear.ran &&
                    filesystem::exists(filesystem::path(cfg.out) / "linear_response_summary.txt");
        ok = ok && pass;
        cout << "Test 4 degree precondition + predictor: " << (pass ? "PASS" : "FAIL")
             << " S_degree=" << out.degree.S_degree
             << " S_pred=" << out.linear.S_pred << '\n';
    }
    {
        Config cfg;
        cfg.nx = 8; cfg.ny = 8; cfg.flips = 50; cfg.seed = 1; cfg.radius_seed = 101;
        cfg.target_dist = "lognormal"; cfg.target_delta = 0.25; cfg.target_radius_ratio = 5.0;
        cfg.compare_topology_strategies = true; cfg.degree_pre_steps = 20; cfg.degree_pre_trials = 10;
        cfg.linear_surgery_steps = 20; cfg.linear_surgery_trials_per_step = 10;
        cfg.strategy_out = "strategy_compare_v8_test.csv"; cfg.out = "strategy_compare_v8_test";
        vector<RunOutput> rows = run_strategy_compare(cfg);
        bool file_exists = filesystem::exists(cfg.strategy_out);
        bool pass = rows.size() == 4 && file_exists;
        ok = ok && pass;
        cout << "Test 5 strategy compare small: " << (pass ? "PASS" : "FAIL")
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

        if (cfg.degree_fraction_scan) {
            double t0 = now_sec();
            vector<RunOutput> rows = run_degree_fraction_scan(cfg);
            int predicted = 0, refined = 0;
            for (const RunOutput& r : rows) {
                if (r.linear.linear_predicted_good) predicted++;
                if (r.packing.physical_INC_like_candidate_v6) refined++;
            }
            cout << "degree_fraction_rows " << rows.size()
                 << " degree_fraction_out " << cfg.degree_fraction_out
                 << " predicted_good " << predicted
                 << " refined_physical " << refined
                 << " total_sec " << scientific << setprecision(6) << (now_sec() - t0) << '\n';
            return 0;
        }

        if (cfg.compare_topology_strategies) {
            double t0 = now_sec();
            vector<RunOutput> rows = run_strategy_compare(cfg);
            cout << "strategy_rows " << rows.size()
                 << " strategy_out " << cfg.strategy_out
                 << " total_sec " << scientific << setprecision(6) << (now_sec() - t0) << '\n';
            return 0;
        }

        if (cfg.screen_topologies) {
            double t0 = now_sec();
            run_topology_screen(cfg);
            cout << "screen_rows " << cfg.screen_count
                 << " screen_out " << cfg.screen_out
                 << " total_sec " << scientific << setprecision(6) << (now_sec() - t0) << '\n';
            return 0;
        }

        if (cfg.linear_pareto_scan) {
            double t0 = now_sec();
            vector<RunOutput> rows = run_linear_pareto_scan(cfg);
            int predicted_good = 0, refined_physical = 0;
            for (const RunOutput& r : rows) {
                if (r.linear.linear_predicted_good) predicted_good++;
                if (r.packing.physical_INC_like_candidate_v6) refined_physical++;
            }
            cout << "linear_pareto_rows " << rows.size()
                 << " pareto_out " << cfg.pareto_out
                 << " predicted_good " << predicted_good
                 << " refined_physical " << refined_physical
                 << " total_sec " << scientific << setprecision(6) << (now_sec() - t0) << '\n';
            return 0;
        }

        if (cfg.pareto_scan) {
            double t0 = now_sec();
            vector<RunOutput> rows = run_pareto_scan(cfg);
            double total = now_sec() - t0;
            int fixed_good = 0, valid_disk = 0, strong = 0, physical_v6 = 0;
            for (const RunOutput& r : rows) {
                if (r.packing.fixed_intrinsic_good) fixed_good++;
                if (r.packing.valid_disk_packing_candidate) valid_disk++;
                if (r.packing.strong_INC_candidate) strong++;
                if (r.packing.physical_INC_like_candidate_v6) physical_v6++;
            }
            cout << "pareto_rows " << rows.size()
                 << " pareto_out " << cfg.pareto_out
                 << " fixed_good " << fixed_good
                 << " valid_disk " << valid_disk
                 << " strong " << strong
                 << " physical_v6 " << physical_v6
                 << " total_sec " << scientific << setprecision(6) << total << '\n';
            return 0;
        }

        if (cfg.multistart_surgery) {
            double t0 = now_sec();
            vector<RunOutput> rows = run_multistart_surgery(cfg);
            double total = now_sec() - t0;
            int valid_disk = 0, strong = 0, physical = 0, physical_v6 = 0;
            for (const RunOutput& r : rows) {
                if (r.packing.valid_disk_packing_candidate) valid_disk++;
                if (r.packing.strong_INC_candidate) strong++;
                if (r.packing.physical_INC_like_candidate) physical++;
                if (r.packing.physical_INC_like_candidate_v6) physical_v6++;
            }
            cout << "multistart_rows " << rows.size()
                 << " multistart_out " << cfg.multistart_out
                 << " valid_disk " << valid_disk
                 << " strong " << strong
                 << " physical " << physical
                 << " physical_v6 " << physical_v6
                 << " total_sec " << scientific << setprecision(6) << total << '\n';
            return 0;
        }

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

        RunOutput out;
        if (cfg.surgery_metric == "linear" && cfg.proposal_mode != "random") {
            out = run_linear_surgery(cfg);
        } else {
            out = cfg.surgery ? run_surgery(cfg) : execute_run(cfg, true, cfg.progress_every > 0);
        }
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
        cout << "final_rms_K " << out.ricci.state.rms_K << '\n';
        cout << "final_E_reg " << out.ricci.final_E_reg << '\n';
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
        cout << "face_degenerate_count " << out.packing.face_degenerate_count << '\n';
        cout << "face_negative_signed_count " << out.packing.face_negative_signed_count << '\n';
        cout << "orientation_neighbor_violation_count " << out.packing.orientation_neighbor_violation_count << '\n';
        cout << "degree6_fraction " << out.packing.degree6_fraction << '\n';
        cout << "valid_disk_packing_candidate " << (out.packing.valid_disk_packing_candidate ? "yes" : "no") << '\n';
        cout << "strong_INC_candidate " << (out.packing.strong_INC_candidate ? "yes" : "no") << '\n';
        cout << "physical_INC_like_candidate " << (out.packing.physical_INC_like_candidate ? "yes" : "no") << '\n';
        cout << "fixed_intrinsic_good " << (out.packing.fixed_intrinsic_good ? "yes" : "no") << '\n';
        cout << "weak_distribution_preserved " << (out.packing.weak_distribution_preserved ? "yes" : "no") << '\n';
        cout << "distribution_preserved " << (out.packing.distribution_preserved ? "yes" : "no") << '\n';
        cout << "physical_INC_like_candidate_v6 " << (out.packing.physical_INC_like_candidate_v6 ? "yes" : "no") << '\n';
        cout << "D_sort " << out.dist.D_sort << '\n';
        cout << "D_quantile " << out.dist.D_quantile << '\n';
        cout << "D_delta " << out.dist.D_delta << '\n';
        cout << "D_ratio " << out.dist.D_ratio << '\n';
        if (out.degree.ran) {
            cout << "S_degree " << out.degree.S_degree << '\n';
            cout << "E_deg " << out.degree.E_deg << '\n';
            cout << "max_degree_error " << out.degree.max_degree_error << '\n';
            cout << "degree_match_fraction " << out.degree.degree_match_fraction << '\n';
            cout << "degree_radius_corr " << out.degree.degree_radius_corr << '\n';
            cout << "target_degree_radius_corr " << out.degree.target_degree_radius_corr << '\n';
        }
        cout << "time_topology_sec " << out.time.topology << '\n';
        cout << "time_ricci_sec " << out.time.ricci << '\n';
        cout << "time_develop_sec " << out.time.develop << '\n';
        cout << "time_total_sec " << out.time.total << '\n';
        cout << "stopped_reason " << out.ricci.stopped_reason << '\n';
        if (out.linear.ran) {
            cout << "S_pred " << out.linear.S_pred << '\n';
            cout << "S_total " << linear_total_score(out, cfg) << '\n';
            cout << "rms_K0 " << out.linear.rms_K0 << '\n';
            cout << "max_K0 " << out.linear.max_K0 << '\n';
            cout << "rms_K_lin " << out.linear.rms_K_lin << '\n';
            cout << "max_K_lin " << out.linear.max_K_lin << '\n';
            cout << "du_rms " << out.linear.du_rms << '\n';
            cout << "du_max " << out.linear.du_max << '\n';
            cout << "degree_radius_corr " << out.linear.degree_radius_corr << '\n';
            cout << "singular_response_flag " << (out.linear.singular_response_flag ? "yes" : "no") << '\n';
            cout << "linear_predicted_good " << (out.linear.linear_predicted_good ? "yes" : "no") << '\n';
        }
        return 0;
    } catch (const exception& e) {
        cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
