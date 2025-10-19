#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "gf_utils.h"
#include "graphflood/define_types.h"
#include "topotoolbox.h"
#include "pq_maxheap.h"

/*
 * GRAPHFLOOD FULL IMPLEMENTATION
 *
 * This module implements the complete GraphFlood algorithm for transient flood
 * simulation. Unlike the metrics function which calculates steady-state flow,
 * this implementation solves the time-dependent shallow water equations using
 * an iterative approach.
 *
 * The algorithm supports both:
 * - Single Flow Direction (SFD): Flow follows steepest descent path
 * - Multiple Flow Direction (MFD): Flow is distributed proportionally to all
 * downslope neighbors
 *
 * Key Features:
 * - Transient simulation with explicit time stepping
 * - Priority flooding to handle depressions
 * - Mass conservation through continuity equation
 * - Manning's friction for discharge calculations
 * - Flexible boundary conditions
 */

// ============================================================================
// UTILITY FUNCTIONS: Min/Max operations for floating point values
// ============================================================================

// Helper functions to mimic C++ std::min and std::max
inline GF_FLOAT max_float(GF_FLOAT a, GF_FLOAT b) { return (a > b) ? a : b; }
inline GF_FLOAT min_float(GF_FLOAT a, GF_FLOAT b) { return (a < b) ? a : b; }

// ============================================================================
// SINGLE FLOW DIRECTION IMPLEMENTATION
// ============================================================================

/*
 * Internal function implementing GraphFlood with Single Flow Direction (SFD)
 *
 * In SFD mode:
 * - Each cell directs all its flow to exactly one downstream neighbor
 * - Flow direction is determined by steepest hydraulic gradient
 * - More computationally efficient but less realistic for divergent flow areas
 * - Better suited for channelized flow systems
 *
 * Algorithm:
 * 1. Build flow graph using priority flooding
 * 2. Accumulate flow following single steepest paths
 * 3. Calculate discharge using Manning's equation
 * 4. Update water depths using continuity equation
 * 5. Repeat for specified number of iterations
 */
void _graphflood_full_sfd(
    GF_FLOAT* Z,               // Digital elevation model [input/const]
    GF_FLOAT* hw,              // Water depth array [input/output]
    uint8_t* BCs,              // Boundary conditions [input/const]
    GF_FLOAT* Precipitations,  // Precipitation rates [input/const]
    GF_FLOAT* manning,         // Manning's roughness [input/const]
    GF_UINT* dim,              // Grid dimensions [nx, ny] [input/const]
    GF_FLOAT dt,               // Time step size [input]
    GF_FLOAT dx,               // Grid spacing [input]
    bool SFD,                  // Flow direction flag [input]
    bool D8,                   // Connectivity flag [input]
    GF_UINT N_iterations,      // Number of iterations [input]
    GF_FLOAT step)             // Flooding step size [input]
{
  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION: Create working arrays for SFD algorithm
  // --------------------------------------------------------------------------

  // Hydraulic surface elevation (ground + water)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i)
    Zw[i] = Z[i] + hw[i];  // Water surface = elevation + depth

  // Flow graph data structures for single flow direction
  GF_UINT* Sreceivers =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));  // Single receiver per cell
  GF_FLOAT* distToReceivers =
      (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));  // Distance to receiver
  GF_UINT* Sdonors = (GF_UINT*)malloc(
      sizeof(GF_UINT) * nxy(dim) * (D8 ? 8 : 4));  // Multiple donors per cell
  uint8_t* NSdonors = (uint8_t*)malloc(sizeof(uint8_t) *
                                       nxy(dim));  // Number of donors per cell
  GF_UINT* Stack =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));  // Processing order stack
  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) *
                                     nxy(dim));  // Input discharge per cell

  // Cell area for volume calculations
  GF_FLOAT cell_area = dx * dx;

  // --------------------------------------------------------------------------
  // MAIN ITERATION LOOP: Time stepping for transient simulation
  // --------------------------------------------------------------------------

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // ------------------------------------------------------------------------
    // STEP 1: BUILD FLOW GRAPH with priority flooding
    // ------------------------------------------------------------------------
    /*
     * Priority flooding ensures:
     * - All depressions in hydraulic surface are filled
     * - Flow paths are monotonically decreasing (no local minima)
     * - Single receiver is identified for each cell
     * - Topological ordering is established
     *
     * This step is crucial for numerical stability and physical realism
     */
    compute_sfgraph_priority_flood(Zw, Sreceivers, distToReceivers, Sdonors,
                                   NSdonors, Stack, BCs, dim, dx, D8, step);

    // ------------------------------------------------------------------------
    // STEP 2: FLOW ACCUMULATION following single flow paths
    // ------------------------------------------------------------------------
    /*
     * Accumulate discharge from precipitation and upstream contributions
     * Following single steepest descent paths established by flow graph
     */
    compute_weighted_drainage_area_single_flow(Qwin, Precipitations, Sreceivers,
                                               Stack, dim, dx);

    // ------------------------------------------------------------------------
    // STEP 3: UPDATE WATER DEPTHS using continuity equation
    // ------------------------------------------------------------------------
    /*
     * Process cells in reverse topological order (upstream to downstream)
     * This ensures hydraulic gradients remain explicit during updates
     */
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      // Get current cell (reverse stack order for proper upstream-downstream
      // processing)
      GF_UINT node = Stack[nxy(dim) - i - 1];
      GF_UINT rec = Sreceivers[node];

      // Skip cells that flow to themselves (pits/outlets)
      if (rec == node) continue;

      // Skip boundary cells that can discharge out of domain
      if (can_out(node, BCs)) continue;

      // Skip dry cells with no input (computational efficiency)
      if (Zw[node] == Z[node] && Qwin[node] == 0) continue;

      // ----------------------------------------------------------------------
      // HYDRAULIC GRADIENT CALCULATION
      // ----------------------------------------------------------------------

      // Calculate water surface slope to receiver
      // Use minimum slope to ensure numerical stability
      GF_FLOAT tSw =
          min_float(Zw[node] - Zw[rec], (GF_FLOAT)1e-6) / distToReceivers[node];

      // ----------------------------------------------------------------------
      // MANNING'S DISCHARGE CALCULATION
      // ----------------------------------------------------------------------

      /*
       * Manning's equation for open channel flow:
       * Q = (1/n) × A × R^(2/3) × S^(1/2)
       *
       * For 2D overland flow approximation:
       * Q ≈ (width/n) × depth^(5/3) × slope^(1/2)
       */
      GF_FLOAT tQwout = 0.0;
      if (Zw[node] > Z[node]) {  // Only if water is present
        GF_FLOAT depth = Zw[node] - Z[node];
        tQwout = (GF_FLOAT)(distToReceivers[node] / manning[node] *
                            pow(depth, 5.0 / 3.0) * sqrt(tSw));
      }

      // ----------------------------------------------------------------------
      // CONTINUITY EQUATION: Update water surface elevation
      // ----------------------------------------------------------------------

      /*
       * Explicit finite difference for continuity equation:
       * ∂h/∂t = (Q_in - Q_out) / A_cell
       *
       * h^(n+1) = h^n + dt × (Q_in - Q_out) / A_cell
       *
       * Ensures water surface never goes below ground level
       */
      Zw[node] =
          max_float(Z[node], Zw[node] + dt * (Qwin[node] - tQwout) / cell_area);
    }
  }

  // --------------------------------------------------------------------------
  // FINALIZATION: Convert hydraulic surface back to water depth
  // --------------------------------------------------------------------------

  // Extract water depths from hydraulic surface
  for (GF_UINT i = 0; i < nxy(dim); ++i)
    hw[i] = max_float(0.0, Zw[i] - Z[i]);  // Ensure non-negative depths

  // Free allocated memory
  free(Zw);
  free(Qwin);
  free(Sreceivers);
  free(distToReceivers);
  free(Sdonors);
  free(NSdonors);
  free(Stack);
}

// ============================================================================
// MULTIPLE FLOW DIRECTION IMPLEMENTATION
// ============================================================================

/*
 * Internal function implementing GraphFlood with Multiple Flow Direction (MFD)
 *
 * In MFD mode:
 * - Each cell can direct flow to multiple downstream neighbors
 * - Flow is distributed proportionally based on hydraulic gradients
 * - More physically realistic for divergent flow areas (alluvial fans,
 * floodplains)
 * - Computationally more intensive due to flow splitting calculations
 *
 * Algorithm:
 * 1. Apply priority flooding to hydraulic surface
 * 2. For each cell, calculate gradients to all valid neighbors
 * 3. Distribute flow proportionally based on gradient weights
 * 4. Calculate discharge using Manning's equation along steepest path
 * 5. Update water depths using continuity equation
 * 6. Repeat for specified iterations
 */
void _graphflood_full_mfd(
    GF_FLOAT* Z,               // Digital elevation model [input/const]
    GF_FLOAT* hw,              // Water depth array [input/output]
    uint8_t* BCs,              // Boundary conditions [input/const]
    GF_FLOAT* Precipitations,  // Precipitation rates [input/const]
    GF_FLOAT* manning,         // Manning's roughness [input/const]
    GF_UINT* dim,              // Grid dimensions [nx, ny] [input/const]
    GF_FLOAT dt,               // Time step size [input]
    GF_FLOAT dx,               // Grid spacing [input]
    bool SFD,                  // Flow direction flag [input]
    bool D8,                   // Connectivity flag [input]
    GF_UINT N_iterations,      // Number of iterations [input]
    GF_FLOAT step)             // Flooding step size [input]
{
  // --------------------------------------------------------------------------
  // NEIGHBOR CONNECTIVITY SETUP
  // --------------------------------------------------------------------------

  // Initialize neighbor offset arrays for grid traversal
  GF_INT offset[8];
  (D8 == false)
      ? generate_offset_D4_flat(offset, dim)  // 4-connectivity (cardinal)
      : generate_offset_D8_flat(offset,
                                dim);  // 8-connectivity (cardinal + diagonal)

  // Initialize distance arrays for each neighbor direction
  GF_FLOAT offdx[8];
  (D8 == false)
      ? generate_offsetdx_D4(offdx, dx)   // Cardinal distances only
      : generate_offsetdx_D8(offdx, dx);  // Cardinal + diagonal distances

  // Diagonal distance and cell area for calculations
  GF_FLOAT dxy = (GF_FLOAT)sqrt(2) * dx;  // Diagonal distance
  GF_FLOAT cell_area = dx * dx;           // Cell area for volume calculations

  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION: Create working arrays for MFD algorithm
  // --------------------------------------------------------------------------

  // Hydraulic surface elevation (ground + water)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i) Zw[i] = Z[i] + hw[i];

  GF_FLOAT* Qwin =
      (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));  // Input discharge
  GF_FLOAT* Qwout =
      (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));  // Output discharge
  GF_UINT* Stack =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));  // Processing order

  // Initialize arrays
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qwin[i] = 0.0;   // Start with zero input discharge
    Qwout[i] = 0.0;  // Start with zero output discharge
    Stack[i] = i;    // Initial stack order
  }

  // Flow weight array for multiple flow direction calculations
  GF_FLOAT weights[8];

  // --------------------------------------------------------------------------
  // MAIN ITERATION LOOP: Time stepping for transient simulation
  // --------------------------------------------------------------------------

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // ------------------------------------------------------------------------
    // STEP 1: PRIORITY FLOODING and topological ordering
    // ------------------------------------------------------------------------
    /*
     * Priority flooding fills depressions and establishes flow paths
     * Topological ordering ensures proper upstream-downstream processing
     */
    compute_priority_flood_plus_topological_ordering(Zw, Stack, BCs, dim, D8,
                                                     step);

    // Reset discharge arrays for this iteration
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      Qwin[i] = 0.0;
      Qwout[i] = 0.0;
    }

    // ------------------------------------------------------------------------
    // STEP 2: FLOW ACCUMULATION with multiple flow direction
    // ------------------------------------------------------------------------
    /*
     * Process cells in reverse topological order (upstream to downstream)
     * This maintains explicit hydraulic gradients during flow distribution
     */
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      // Get current cell index (reverse order processing)
      GF_UINT node = Stack[nxy(dim) - i - 1];

      // Skip invalid cells
      if (is_nodata(node, BCs)) continue;

      // Skip boundary cells that discharge out of domain
      if (can_out(node, BCs)) continue;

      // ----------------------------------------------------------------------
      // LOCAL PRECIPITATION INPUT
      // ----------------------------------------------------------------------

      // Add precipitation contribution (rate × cell area = volume/time)
      Qwin[node] += Precipitations[node] * dx * dx;

      // ----------------------------------------------------------------------
      // GRADIENT ANALYSIS: Calculate slopes to all valid neighbors
      // ----------------------------------------------------------------------

      GF_FLOAT sumslope = 0.0;  // Total weighted slope for flow partitioning
      GF_FLOAT maxslope = 0.0;  // Steepest slope for Manning's calculation
      GF_FLOAT dxmaxdir = dx;   // Distance to steepest neighbor

      // Examine all potential neighbors
      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        // Check if neighbor is valid and within bounds
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
          weights[n] = 0;
          continue;
        }

        GF_UINT nnode = node + offset[n];  // Neighbor cell index

        // Skip uphill neighbors or flow-restricted cells
        if (Zw[nnode] >= Zw[node] ||             // Neighbor is uphill
            can_receive(nnode, BCs) == false ||  // Neighbor cannot receive flow
            can_give(node, BCs) == false) {  // Current cell cannot give flow
          weights[n] = 0;
          continue;
        }

        // Calculate hydraulic gradient (slope) to this neighbor
        // Ensure minimum slope to prevent numerical issues
        GF_FLOAT tSw =
            max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);

        // Calculate flow weight proportional to gradient and effective flow
        // width Use appropriate distance based on flow direction (cardinal vs
        // diagonal)
        weights[n] = tSw * ((dx == offdx[n] || D8 == false) ? dx : dxy);

        // Accumulate total slope for proportional distribution
        sumslope += weights[n];

        // Track steepest gradient for discharge calculations
        if (tSw > maxslope) {
          maxslope = tSw;
          dxmaxdir = offdx[n];
        }
      }

      // ----------------------------------------------------------------------
      // FLOW DISTRIBUTION: Partition flow among downstream neighbors
      // ----------------------------------------------------------------------

      /*
       * Multiple Flow Direction (MFD) algorithm:
       * Each downstream neighbor receives flow proportional to its gradient
       * weight Flow fraction = (neighbor_weight / total_weights)
       *
       * This approach naturally handles:
       * - Flow convergence (multiple inputs to one cell)
       * - Flow divergence (one cell distributing to multiple neighbors)
       * - Flow routing through complex topography
       */
      if (sumslope > 0) {
        for (GF_UINT n = 0; n < N_neighbour(D8); ++n) {
          if (weights[n] == 0) continue;  // Skip neighbors with no flow

          // Distribute flow proportionally to gradient weight
          Qwin[node + offset[n]] += weights[n] / sumslope * Qwin[node];
        }
      }

      // ----------------------------------------------------------------------
      // MANNING'S DISCHARGE CALCULATION
      // ----------------------------------------------------------------------

      /*
       * Calculate discharge using Manning's equation along steepest gradient
       * This represents the maximum conveyance capacity of the cell
       */
      if (Zw[node] > Z[node]) {  // Only if water is present
        GF_FLOAT depth = Zw[node] - Z[node];
        Qwout[node] = (GF_FLOAT)(dxmaxdir / manning[node] *
                                 pow(depth, 5.0 / 3.0) * sqrt(maxslope));
      }
    }

    // ------------------------------------------------------------------------
    // STEP 3: UPDATE WATER DEPTHS using continuity equation
    // ------------------------------------------------------------------------

    for (GF_UINT node = 0; node < nxy(dim); ++node) {
      // ----------------------------------------------------------------------
      // CONTINUITY EQUATION: Update hydraulic surface
      // ----------------------------------------------------------------------

      /*
       * Explicit time stepping for continuity equation:
       *
       * ∂(h)/∂t + ∇·(q) = P
       *
       * Where:
       * - h = water depth
       * - q = discharge per unit width
       * - P = precipitation rate
       *
       * Discretized as:
       * h^(n+1) = h^n + dt × (Q_in - Q_out) / A_cell
       *
       * Constraint: Water surface cannot go below ground level
       */
      Zw[node] = max_float(
          Z[node], Zw[node] + dt * (Qwin[node] - Qwout[node]) / cell_area);
    }
  }

  // --------------------------------------------------------------------------
  // FINALIZATION: Convert hydraulic surface back to water depths
  // --------------------------------------------------------------------------

  // Extract water depths, ensuring non-negative values
  for (GF_UINT i = 0; i < nxy(dim); ++i) hw[i] = max_float(0.0, Zw[i] - Z[i]);

  // Free allocated memory
  free(Zw);
  free(Qwin);
  free(Qwout);
  free(Stack);
}

// ============================================================================
// PUBLIC API FUNCTION: Main GraphFlood interface
// ============================================================================

/*
 * GRAPHFLOOD_FULL: Complete transient flood simulation
 *
 * This is the main public interface for running complete flood simulations.
 * It supports both steady-state and transient modeling with flexible
 * flow direction schemes and boundary conditions.
 *
 * PARAMETERS:
 * -----------
 * Z            : Digital elevation model [m] (input, constant)
 * hw           : Water depth at each cell [m] (input/output, modified)
 * BCs          : Boundary condition flags (input, constant)
 * Precipitations: Precipitation rate per cell [m/s] (input, constant)
 * manning      : Manning's roughness coefficient [s/m^(1/3)] (input, constant)
 * dim          : Grid dimensions [nx, ny] (input, constant)
 * dt           : Time step size [s] (input, constant)
 * dx           : Grid cell spacing [m] (input, constant)
 * SFD          : Single Flow Direction flag (true=SFD, false=MFD) (input)
 * D8           : Connectivity scheme (true=8-neighbor, false=4-neighbor)
 * (input) N_iterations : Number of time steps to simulate (input) step : Step
 * size for priority flooding algorithm (input)
 *
 * ALGORITHM SELECTION:
 * -------------------
 * The function automatically selects the appropriate algorithm based on the SFD
 * flag:
 *
 * SFD=true:  Single Flow Direction
 *            - Faster computation
 *            - Good for channelized flow
 *            - Each cell flows to exactly one neighbor
 *
 * SFD=false: Multiple Flow Direction
 *            - More realistic for overland flow
 *            - Better for complex topography
 *            - Flow can be distributed among multiple neighbors
 *
 * NUMERICAL CONSIDERATIONS:
 * ------------------------
 * - Time step (dt) should satisfy CFL condition for numerical stability
 * - Grid resolution (dx) affects computational cost and accuracy
 * - Number of iterations determines simulation duration
 * - Priority flooding ensures numerical stability in depressions
 *
 * APPLICATIONS:
 * ------------
 * - Urban flood modeling
 * - Dam break simulations
 * - Overland flow routing
 * - Watershed hydrology
 * - Emergency response planning
 */
TOPOTOOLBOX_API
void graphflood_full(GF_FLOAT* Z,   // Digital elevation model [input]
                     GF_FLOAT* hw,  // Water depths [input/output]
                     uint8_t* BCs,  // Boundary conditions [input]
                     GF_FLOAT* Precipitations,  // Precipitation rates [input]
                     GF_FLOAT* manning,         // Manning's roughness [input]
                     GF_UINT* dim,              // Grid dimensions [input]
                     GF_FLOAT dt,               // Time step size [input]
                     GF_FLOAT dx,               // Grid spacing [input]
                     bool SFD,                  // Flow direction scheme [input]
                     bool D8,                   // Connectivity scheme [input]
                     GF_UINT N_iterations,      // Number of iterations [input]
                     GF_FLOAT step)             // Flooding step size [input]
{
  // Route to appropriate algorithm based on flow direction scheme
  if (SFD) {
    // Single Flow Direction: computationally efficient
    _graphflood_full_sfd(Z, hw, BCs, Precipitations, manning, dim, dt, dx, SFD,
                         D8, N_iterations, step);
  } else {
    // Multiple Flow Direction: more physically realistic
    _graphflood_full_mfd(Z, hw, BCs, Precipitations, manning, dim, dt, dx, SFD,
                         D8, N_iterations, step);
  }
}

/*
 * GRAPHFLOOD METRICS FUNCTION
 *
 * This function implements a hydrological flow routing algorithm using priority
 * flooding and topological ordering to calculate various discharge metrics for
 * flood modeling.
 *
 * The algorithm performs:
 * 1. Priority flooding to determine hydraulic surfaces and flow directions
 * 2. Topological ordering to ensure upstream-to-downstream processing
 * 3. Flow accumulation based on precipitation inputs
 * 4. Discharge calculations using Manning's equation
 *
 * Key metrics calculated:
 * - Qi: Stationary discharge based on drainage area accumulation
 * - Qo: Calculated discharge from the model (should equal Qi at convergence)
 * - q: Discharge per unit width
 * - u: Flow velocity
 * - Sw: Water surface slope
 */

TOPOTOOLBOX_API
void graphflood_metrics(
    GF_FLOAT* Z,               // Digital elevation model (DEM) [input]
    GF_FLOAT* hw,              // Water depth at each cell [input]
    uint8_t* BCs,              // Boundary conditions array [input]
    GF_FLOAT* Precipitations,  // Precipitation rate per cell [input]
    GF_FLOAT* manning,         // Manning's roughness coefficient [input]
    GF_UINT* dim,              // Grid dimensions [nx, ny] [input]
    GF_FLOAT dx,               // Grid cell size (spacing) [input]
    bool D8,        // Flow direction scheme: true=D8, false=D4 [input]
    GF_FLOAT step,  // Step size for flooding algorithm [input]
    GF_FLOAT* Qi,   // Input discharge (accumulated flow) [output]
    GF_FLOAT* Qo,   // Output discharge (Manning's equation) [output]
    GF_FLOAT* qo,   // Discharge per unit width [output]
    GF_FLOAT* u,    // Flow velocity [output]
    GF_FLOAT* Sw)   // Water surface slope [output]
{
  // ============================================================================
  // INITIALIZATION: Set up neighbor connectivity and distance arrays
  // ============================================================================

  // Initialize offset array for neighboring cell operations
  // D4: 4 cardinal directions (N,S,E,W), D8: 8 directions (including diagonals)
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)   // 4-connectivity
                : generate_offset_D8_flat(offset, dim);  // 8-connectivity

  // Initialize distance to each neighbor
  // Cardinal directions = dx, diagonal directions = dx*sqrt(2)
  GF_FLOAT offdx[8];
  (D8 == false)
      ? generate_offsetdx_D4(offdx, dx)   // Only cardinal distances
      : generate_offsetdx_D8(offdx, dx);  // Cardinal + diagonal distances

  // Diagonal distance for flow calculations
  GF_FLOAT dxy = (GF_FLOAT)sqrt(2) * dx;

  // ============================================================================
  // MEMORY ALLOCATION: Create working arrays
  // ============================================================================

  // Create hydraulic surface array (elevation + water depth)
  // Zw represents the water surface elevation at each cell
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i)
    Zw[i] = Z[i] + hw[i];  // Water surface = ground elevation + water depth

  // Stack for topologically ordered processing of cells
  // Will contain cell indices in processing order (upstream to downstream)
  GF_UINT* Stack = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));

  // Initialize input discharge array to zero and populate stack with cell
  // indices
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qi[i] = 0.0;   // Start with no accumulated flow
    Stack[i] = i;  // Initially, stack contains all cells in order
  }

  // Array to store flow weights for multi-direction flow distribution
  GF_FLOAT weights[8];

  // ============================================================================
  // PRIORITY FLOODING: Determine flow directions and processing order
  // ============================================================================

  /*
   * Priority flooding algorithm:
   * 1. Fills depressions in the hydraulic surface (Zw)
   * 2. Ensures monotonic flow paths (no local minima except outlets)
   * 3. Creates topological ordering in Stack for upstream-to-downstream
   * processing
   *
   * This is critical because it ensures that when we process flow accumulation,
   * all upstream cells are processed before downstream cells.
   */
  compute_priority_flood_plus_topological_ordering(Zw, Stack, BCs, dim, D8,
                                                   step);

  // Re-initialize input discharge array after flooding
  // (flooding may have modified the processing order)
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qi[i] = 0.0;
  }

  // ============================================================================
  // FLOW ACCUMULATION: Process cells in topological order
  // ============================================================================

  /*
   * Process nodes from downstream to upstream (reverse stack order)
   * This ensures that receivers are never processed before donors,
   * maintaining explicit hydraulic slopes throughout the calculation
   */
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    // Get current cell index (processing in reverse topological order)
    GF_UINT node = Stack[nxy(dim) - i - 1];

    // Skip cells with no data
    if (is_nodata(node, BCs)) continue;

    // Skip boundary cells that can flow out of the domain
    // These cells don't accumulate flow from upstream
    if (can_out(node, BCs)) continue;

    // ------------------------------------------------------------------------
    // PRECIPITATION INPUT: Add local water input
    // ------------------------------------------------------------------------

    // Add precipitation contribution to input discharge
    // Precipitation rate × cell area gives volumetric flow rate
    Qi[node] += Precipitations[node] * dx * dx;

    // ------------------------------------------------------------------------
    // FLOW DIRECTION ANALYSIS: Calculate gradients to all neighbors
    // ------------------------------------------------------------------------

    GF_FLOAT sumslope = 0.0;  // Sum of all downslope gradients (for weighting)
    GF_FLOAT maxslope = 0.0;  // Steepest gradient (for Manning's equation)
    GF_FLOAT dxmaxdir = dx;   // Distance to steepest neighbor

    // Examine all possible neighbors
    for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
      // Check if neighbor is within grid bounds and valid
      if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
        weights[n] = 0;  // No flow to invalid neighbors
        continue;
      }

      GF_UINT nnode = node + offset[n];  // Calculate neighbor index

      // Skip uphill neighbors or neighbors that can't receive/give flow
      if (Zw[nnode] >= Zw[node] ||             // Uphill neighbor
          can_receive(nnode, BCs) == false ||  // Neighbor can't receive flow
          can_give(node, BCs) == false) {      // Current cell can't give flow
        weights[n] = 0;
        continue;
      }

      // Calculate hydraulic gradient (slope) to this neighbor
      // Ensure minimum slope to avoid division by zero
      GF_FLOAT tSw =
          max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);

      // Calculate flow weight proportional to gradient and effective width
      // For diagonal flow, use diagonal distance; for cardinal, use dx
      weights[n] = tSw * ((dx == offdx[n] || D8 == false) ? dx : dxy);

      // Accumulate total slope for proportional flow distribution
      sumslope += weights[n];

      // Track steepest gradient for velocity calculations
      if (tSw > maxslope) {
        maxslope = tSw;
        dxmaxdir = offdx[n];  // Distance to steepest neighbor
      }
    }

    // ------------------------------------------------------------------------
    // FLOW DISTRIBUTION: Transfer accumulated flow to downstream neighbors
    // ------------------------------------------------------------------------

    /*
     * Distribute the accumulated flow (Qi) to downstream neighbors
     * proportionally based on their hydraulic gradients.
     * This implements multiple flow direction (MFD) routing.
     */
    if (sumslope > 0) {
      for (GF_UINT n = 0; n < N_neighbour(D8); ++n) {
        if (weights[n] == 0) continue;  // Skip neighbors with no flow

        // Transfer flow proportional to gradient weight
        // Each downstream neighbor receives: (its_weight / total_weight) ×
        // total_flow
        Qi[node + offset[n]] += weights[n] / sumslope * Qi[node];
      }
    }

    // ------------------------------------------------------------------------
    // DISCHARGE CALCULATIONS: Apply Manning's equation for flow velocity
    // ------------------------------------------------------------------------

    /*
     * Calculate volumetric discharge using Manning's friction equation:
     * Q = (1/n) × A × R^(2/3) × S^(1/2)
     *
     * Where:
     * - n = Manning's roughness coefficient
     * - A = Cross-sectional area (approximated as width × depth)
     * - R = Hydraulic radius (approximated as depth for wide channels)
     * - S = Hydraulic slope (water surface gradient)
     *
     * Simplified for 2D flow: Q = (width/n) × depth^(5/3) × slope^(1/2)
     */
    GF_FLOAT tQwout = 0.0;

    // Only calculate discharge if there's water depth
    if (Zw[node] > Z[node]) {
      GF_FLOAT depth = Zw[node] - Z[node];  // Water depth

      // Manning's equation: Q = (width/n) × depth^(5/3) × slope^(1/2)
      tQwout = (GF_FLOAT)(dxmaxdir / manning[node] * pow(depth, 5.0 / 3.0) *
                          sqrt(maxslope));
    }

    // ------------------------------------------------------------------------
    // OUTPUT METRICS: Store calculated values
    // ------------------------------------------------------------------------

    Qo[node] = tQwout;             // Volumetric discharge [m³/s]
    qo[node] = tQwout / dxmaxdir;  // Discharge per unit width [m²/s]

    // Calculate flow velocity: v = q/h (discharge per unit width ÷ depth)
    if (Zw[node] - Z[node] > 0) {
      u[node] = qo[node] / (Zw[node] - Z[node]);  // Flow velocity [m/s]
    } else {
      u[node] = 0.0;  // No velocity if no water depth
    }

    Sw[node] = maxslope;  // Water surface slope [-]
  }

  // Free allocated memory
  free(Zw);
  free(Stack);
}

// ============================================================================
// DYNAMIC INDUCED GRAPH IMPLEMENTATION
// ============================================================================


/*
 * GRAPHFLOOD_DYNAMIC_GRAPH: Dynamic induced graph flood simulation
 *
 * This function implements a wavefront-based flood routing approach that
 * processes cells in descending elevation order, starting from specified
 * input discharge locations. Unlike graphflood_full which processes the
 * entire grid, this method dynamically builds a flow graph by propagating
 * from high to low elevations.
 *
 * Key Differences from graphflood_full:
 * - Uses max-heap to process in descending hydraulic elevation order
 * - Starts from cells with input_Qw > 0
 * - Dynamically expands processing to downstream neighbors
 * - Upstream flow accumulation calculated on-the-fly
 *
 * Algorithm:
 * 1. Initialize from input_Qw cells
 * 2. Process cells in descending Zw order using max-heap
 * 3. For each cell:
 *    - Check for lower neighbors (raise if needed)
 *    - Accumulate flow from upstream neighbors
 *    - Calculate output discharge to downstream neighbors
 *    - Add unvisited downstream neighbors to queue
 * 4. Update water depths after processing all cells
 *
 * PARAMETERS:
 * -----------
 * Z            : Digital elevation model [m] (input, constant)
 * hw           : Water depth at each cell [m] (input/output, modified)
 * BCs          : Boundary condition flags (input, constant)
 * Precipitations: Precipitation rate per cell [m/s] (input, constant)
 * manning      : Manning's roughness coefficient [s/m^(1/3)] (input, constant)
 * input_Qw     : Input discharge locations/values [m³/s] (input, constant)
 * dim          : Grid dimensions [nx, ny] (input, constant)
 * dt           : Time step size [s] (input, constant)
 * dx           : Grid cell spacing [m] (input, constant)
 * D8           : Connectivity scheme (true=8-neighbor, false=4-neighbor)
 * N_iterations : Number of time steps to simulate (input)
 */
TOPOTOOLBOX_API
void graphflood_dynamic_graph(
    GF_FLOAT* Z,               // Digital elevation model [input]
    GF_FLOAT* hw,              // Water depths [input/output]
    uint8_t* BCs,              // Boundary conditions [input]
    GF_FLOAT* Precipitations,  // Precipitation rates [input]
    GF_FLOAT* manning,         // Manning's roughness [input]
    GF_FLOAT* input_Qw,        // Input discharge locations [input]
    GF_UINT* dim,              // Grid dimensions [input]
    GF_FLOAT dt,               // Time step size [input]
    GF_FLOAT dx,               // Grid spacing [input]
    bool D8,                   // Connectivity scheme [input]
    GF_UINT N_iterations)      // Number of iterations [input]
{
  // --------------------------------------------------------------------------
  // NEIGHBOR CONNECTIVITY SETUP
  // --------------------------------------------------------------------------


  // printf("REACHES HERE0\n");
  // fflush(stdout);

  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);

  GF_FLOAT offdx[8];
  (D8 == false) ? generate_offsetdx_D4(offdx, dx)
                : generate_offsetdx_D8(offdx, dx);

  GF_FLOAT dxy = (GF_FLOAT)sqrt(2) * dx;
  GF_FLOAT cell_area = dx * dx;

  // --------------------------------------------------------------------------
  // MEMORY ALLOCATION
  // --------------------------------------------------------------------------

  GF_UINT tnxy = nxy(dim);

  // Hydraulic surface (elevation + water depth)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  for (GF_UINT i = 0; i < tnxy; ++i) {
    Zw[i] = Z[i] + hw[i];
  }

  // Flow arrays
  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);
  GF_FLOAT* Qwout = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * tnxy);

  // State tracking arrays
  uint8_t* visited = (uint8_t*)malloc(sizeof(uint8_t) * tnxy);
  uint8_t* inPQ = (uint8_t*)malloc(sizeof(uint8_t) * tnxy);

  // Weight array for flow distribution
  GF_FLOAT weights[8];

  // --------------------------------------------------------------------------
  // IDENTIFY INPUT CELLS
  // --------------------------------------------------------------------------

  // printf("REACHES HERE1\n");
  fflush(stdout);

  // Count cells with input discharge
  GF_UINT n_input_cells = 0;
  for (GF_UINT i = 0; i < tnxy; ++i) {
    if (input_Qw[i] > 0.0) {
      n_input_cells++;
    }
  }

  // Create array of input cell indices
  GF_UINT* input_indices =
      (GF_UINT*)malloc(sizeof(GF_UINT) * (n_input_cells + 1));
  GF_UINT idx = 0;
  for (GF_UINT i = 0; i < tnxy; ++i) {
    if (input_Qw[i] > 0.0) {
      input_indices[idx++] = i;
    }
  }
  
  // printf("REACHES HERE2\n");
  // fflush(stdout);

  // --------------------------------------------------------------------------
  // MAIN ITERATION LOOP
  // --------------------------------------------------------------------------

  MaxHeapPQueue pq;
  maxheap_init(&pq, tnxy);

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // Reset arrays
    for (GF_UINT i = 0; i < tnxy; ++i) {
      Qwin[i] = input_Qw[i];
      Qwout[i] = 0.0;
      visited[i] = false;
      inPQ[i] = false;
    }


    // printf("REACHES HERE3\n");
    // fflush(stdout);

    // ------------------------------------------------------------------------
    // INITIALIZE PRIORITY QUEUE with input cells
    // ------------------------------------------------------------------------

    for (GF_UINT i = 0; i < n_input_cells; ++i) {
      GF_UINT node = input_indices[i];
      maxheap_push(&pq, node, Zw[node]);
      inPQ[node] = true;
    }

    // printf("REACHES HERE4\n");
    // fflush(stdout);

    // ------------------------------------------------------------------------
    // PROCESS CELLS IN DESCENDING ELEVATION ORDER
    // ------------------------------------------------------------------------

    while (maxheap_empty(&pq) == false) {
      
      // printf("YOLO\n");
      // fflush(stdout);
      GF_UINT node = maxheap_pop_and_get_key(&pq);
      

      // printf("%zu", (size_t)node);
      // printf("\n");
      // fflush(stdout);

      inPQ[node] = false;
      visited[node] = true;

      // Skip invalid cells
      if (is_nodata(node, BCs)) continue;

      // --------------------------------------------------------------------
      // CHECK FOR LOWER NEIGHBORS and raise if needed
      // --------------------------------------------------------------------

      GF_FLOAT min_neighbor_zw = INFINITY;
      bool has_lower_neighbor = false;
      bool has_any_neighbor = false;

      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

        GF_UINT nnode = node + offset[n];
        if (is_nodata(nnode, BCs)) continue;

        has_any_neighbor = true;

        if (Zw[nnode] < Zw[node]) {
          has_lower_neighbor = true;
        }
        if (Zw[nnode] < min_neighbor_zw) {
          min_neighbor_zw = Zw[nnode];
        }
      }

      // Raise cell if no lower neighbors
      if (has_any_neighbor && has_lower_neighbor == false) {
        Zw[node] = min_neighbor_zw + 1e-3;
      }

      // --------------------------------------------------------------------
      // ACCUMULATE FLOW FROM UPSTREAM NEIGHBORS
      // --------------------------------------------------------------------

      // Add local precipitation
      Qwin[node] += Precipitations[node] * cell_area;


      GF_FLOAT sum_slopes_j = 0.0;
      GF_FLOAT maxslope = 0.0;
      GF_FLOAT dxmaxslope = 0.0;
      // Calculate contribution from upstream neighbors
      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

        GF_UINT nnode = node + offset[n];  // Potential upstream neighbor
        if (is_nodata(nnode, BCs)) continue;

        // Check if this is an upstream neighbor (higher elevation)
        if (Zw[nnode] >= Zw[node]) continue;

        // Calculate slopes from upstream neighbor to all its downstream
        // neighbors
        GF_FLOAT slope_j =
              (Zw[node] - Zw[nnode]) / offdx[n];
        sum_slopes_j += slope_j;
        
      }

      if(sum_slopes_j > 0){
        // Calculate contribution from upstream neighbors
        for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
          if (check_bound_neighbour(node, n, dim, BCs, D8) == false) continue;

          GF_UINT nnode = node + offset[n];  // Potential upstream neighbor
          if (is_nodata(nnode, BCs)) continue;

          // Check if this is an upstream neighbor (higher elevation)
          if (Zw[nnode] >= Zw[node]) continue;

          // Calculate slopes from upstream neighbor to all its downstream
          // neighbors
          GF_FLOAT slope_j =
                max_float((GF_FLOAT)1e-8, (Zw[node] - Zw[nnode]) / offdx[n]);

          if (slope_j>maxslope){
            maxslope = slope_j;
            dxmaxslope = offdx[n];
          }


          Qwin[nnode] += (slope_j / sum_slopes_j) * Qwin[node];

          // Add downstream neighbor to queue if not already there
          if (inPQ[nnode] == false && visited[nnode] == false) {
            maxheap_push(&pq, nnode, Zw[nnode]);
            inPQ[nnode] = true;
          }
          
        }

      // Distribute flow proportionally
        
      }
      else{
        maxheap_push(&pq, node, Zw[node]);
        inPQ[node] = true;
        Zw[node] += 1e-3;
        continue;
      }

      // --------------------------------------------------------------------
      // CALCULATE OUTPUT DISCHARGE TO DOWNSTREAM NEIGHBORS
      // --------------------------------------------------------------------

      // Calculate discharge using Manning's equation
      if (Zw[node] > Z[node]) {
        GF_FLOAT depth = maxslope(GF_FLOAT(Zw[node] - Z[node]),0.);
        Qwout[node] =
            (GF_FLOAT)(dxmaxslope / manning[node] * pow(depth, 5.0 / 3.0) *
                       sqrt(maxslope));
      }
    }

    // ------------------------------------------------------------------------
    // UPDATE WATER DEPTHS
    // ------------------------------------------------------------------------

    for (GF_UINT node = 0; node < tnxy; ++node) {
      // Only update cells with water or discharge
      if (Zw[node] > Z[node] || visited[node] || Qwin[node]>0.) {
        Zw[node] = max_float(
            Z[node], Zw[node] + dt * (Qwin[node] - Qwout[node]) / cell_area);
      }
    }
  }

  // --------------------------------------------------------------------------
  // FINALIZATION
  // --------------------------------------------------------------------------

  // Back-calculate water depths
  for (GF_UINT i = 0; i < tnxy; ++i) {
    hw[i] = max_float(0.0, Zw[i] - Z[i]);
  }

  // Free memory
  maxheap_free(&pq);
  free(Zw);
  free(Qwin);
  free(Qwout);
  free(visited);
  free(inPQ);
  free(input_indices);
}
