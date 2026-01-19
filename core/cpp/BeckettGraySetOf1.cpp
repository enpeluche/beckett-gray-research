using namespace std;

#include <fstream>  // File stream: used for reading from or writing to files on disk
#include <iostream> // Input/Output stream: standard console printing (cout, cerr)
#include <vector>   // Dynamic array container: provides resizable storage for data
#include <cmath>    // Mathematical functions: provides pow, sqrt, log, etc.
#include <cstdio>   // C-style Input/Output: used for fast I/O like printf or scanf
#include <cstdlib>  // Standard library: general purpose functions (memory, rand, exit)
#include <ctime>    // C-style Time: used for basic timing and seeding random numbers
#include <bitset>   // Bit manipulation: useful for fixed-size binary representations

#include "EvalPerf.hpp"

// --- Core Algorithm Variables ---
int n;       // Number of bits (dimension of the hypercube)
int n2;      // Total number of nodes (2^n)
int s;       // Current sequence length during recursion
int smax;    // Maximum sequence length reached so far (tracking progress)
int nbAppel; // Recursive function call counter (performance monitoring)

// --- Sequence Storage ---
vector<int> bgc;                  // Current Beckett Gray Code sequence being constructed
vector<vector<int>> listeCodeBGC; // Collection of all valid Beckett Gray Codes found

// --- State Tracking & Constraints ---
vector<bool> avail; // Visitation map (True if the binary integer is unvisited)
vector<int> old;    // FIFO-like buffer of bit indices currently set to 1
vector<int> oldS;   // Timestamps (step indices) when bits were flipped to 1
vector<int> nbMax;  // Balance constraints: max allowed consecutive 1s per bit position

// --- Precomputed Lookup Tables (Optimization) ---
vector<vector<int>> bi;    // Precomputed binary representation for each integer < 2^n
vector<vector<int>> flip;  // Adjacency matrix: flip[i][j] is the neighbor of i by toggling bit j
vector<int> nombreVoisins; // Current degree of each vertex (number of unvisited neighbors)
vector<int> V;             // Hamming weight: V[i] stores the number of set bits (1s) in i
vector<int> nbEltEns;      // Population count: number of remaining elements for each Hamming weight
vector<int> listeVoisinV;  // Edge count between Hamming layers (edges connecting weight k to k+1)

// --- Performance & Statistics ---
EvalPerf PE; // Performance evaluation tool (timers and hardware counters)
int bmin;    // Minimum "max-balance" found across all generated codes
int bmax;    // Maximum "max-balance" found across all generated codes

/**
 * Converts an integer into a binary representation as a vector of bits.
 * @param x The integer to convert.
 * @param n The number of bits to extract (starting from the right).
 * @return A vector of n integers (0 or 1).
 */
vector<int> Binary(int x, int n)
{
    // Convert x to a 32-bit binary string representation
    string b = bitset<32>(static_cast<unsigned long long>(x)).to_string();

    vector<int> res;
    // Extract the last n bits by iterating from the (32-n)th character
    for (int i = b.size() - n; i < b.size(); i++)
    {
        // Convert char ('0' or '1') to int by subtracting the ASCII value of '0'
        res.push_back(b[i] - '0');
    }
    return res;
}

/**
 * Populates the Hamming weight vector (V) and counts elements per weight.
 * @param n The number of bits (dimension).
 * @return A vector where res[i] is the number of set bits (1s) in integer i.
 */
vector<int> listeDes1(int n)
{
    vector<int> res;

    for (int i = 0; i < pow(2, n); i++)
    {
        int cmpt = 0;
        int nb = i;

        // Brian Kernighan's-like bit counting algorithm
        while (nb > 0)
        {
            cmpt += nb & 1; // Check if the least significant bit is 1
            nb >>= 1;       // Right shift to process the next bit
        }
        res.push_back(cmpt);

        // Increment the population count for this specific Hamming weight
        nbEltEns[cmpt]++;
    }
    return res;
}

/**
 * Recursively computes the number of edges between Hamming layers.
 * * In an n-dimensional hypercube, this represents the connections
 * between nodes with weight 'k' and nodes with weight 'k+1'.
 * * @param n The dimension of the hypercube.
 * @return A vector containing the counts of edges between successive layers.
 */
vector<int> pascal(int n)
{
    // Base case for a 1D hypercube (2 nodes, 1 edge connecting them)
    if (n == 1)
    {
        return {2};
    }
    else
    {
        // Recursive step: compute values for n-1
        vector<int> pm1 = pascal(n - 1);
        vector<int> res = {2};

        // Compute internal values using the Pascal triangle addition rule
        for (int i = 1; i < n - 1; ++i)
        {
            res.push_back(pm1[i - 1] + pm1[i]);
        }

        res.push_back(2);
        return res;
    }
}

/**
 * Generates a lookup table containing binary representations of all integers up to 2^n - 1.
 * @param n The number of bits for each representation.
 * @return A 2D vector where res[i] is the binary vector (0s and 1s) of integer i.
 */
vector<vector<int>> listeBinaire(int n)
{
    vector<vector<int>> res;
    for (int i = 0; i < pow(2, n); i++)
    {
        // Compute and store the fixed-length bit vector for each integer
        res.push_back(Binary(i, n));
    }
    return res;
}

/**
 * Checks if a given integer is a power of two using bitwise logic.
 * @param n The integer to check.
 * @return True if n is a power of two, False otherwise.
 */
bool pw2(int n)
{
    // A power of two has exactly one bit set to 1.
    // The operation (n & (n - 1)) clears the lowest set bit.
    // If the result is 0 and n was not 0, then n is a power of two.
    return n != 0 && (n & (n - 1)) == 0;
}

/**
 * Generates an adjacency list for all integers in the 2^n space.
 * * For each integer i, it computes all neighbors reachable by flipping exactly one bit.
 * This effectively maps the edges of an n-dimensional hypercube.
 *
 * @param n The number of bits (dimension of the hypercube).
 * @return A 2D vector where res[i][j] is the result of flipping the j-th bit of i.
 */
vector<vector<int>> listeFlip(int n)
{
    vector<vector<int>> res;

    for (int i = 0; i < pow(2, n); i++)
    {
        res.push_back(vector<int>());
        for (int j = 0; j < n; j++)
        {
            // Bitwise XOR with a single-bit mask (1 << j) toggles the j-th bit
            res[i].push_back(i ^ (1 << j));
        }
    }
    return res;
}

/**
 * Prints a boolean vector to the standard output.
 * * @param l The vector of booleans to be displayed.
 */
void toStringBool(vector<bool> l)
{
    // Iterate through the vector and print 1 for True, 0 for False
    for (bool b : l)
    {
        cout << b << " ";
    }
    cout << endl; // New line after printing all elements
}

/**
 * Prints an integer vector to the standard output.
 * * @param l The vector of integers to be displayed.
 */
void toStringInt(vector<int> l)
{
    // Iterate through the vector and print each integer separated by a space
    for (int i : l)
    {
        cout << i << " ";
    }
    // Print a newline character at the end of the sequence
    cout << endl;
}

/**
 * Analyzes a specific Gray Code to find its maximum run-length of consecutive 1s.
 * Updates the global bmin and bmax if a new best or worst balance is found.
 * * @param i The index of the Gray Code to analyze within listeCodeBGC.
 */
void compteLongeur1(int i)
{
    int max = 0;

    // Check each bit position (column) k
    for (int k = 0; k < n; k++)
    {
        int cmpt = 0;

        // Traverse the entire sequence (n2 = 2^n)
        for (int j = 0; j < n2; j++)
        {
            // Access the k-th bit of the j-th number in the sequence
            if (bi[listeCodeBGC[i][j]][k] == 1)
            {
                cmpt++;
            }
            else
            {
                // Sequence of 1s broken, update the max if necessary
                if (cmpt > max)
                {
                    max = cmpt;
                }
                cmpt = 0;
            }
        }
    }

    // Update global maximum balance if this code has a longer run of 1s
    if (max > bmax)
    {
        cout << i << "max : " << max << endl;
        toStringInt(listeCodeBGC[i]);
        bmax = max;
    };

    // Update global minimum balance if this code is more balanced (shorter max run)
    if (max < bmin)
    {
        cout << i << "min : " << max << endl;
        toStringInt(listeCodeBGC[i]);
        bmin = max;
    }
}

/**
 * Recursive Backtracking function to find Beckett Gray Codes.
 * It explores the hypercube while maintaining the FIFO constraint (Beckett)
 * and pruning the search tree using graph connectivity heuristics.
 * * @param x The current integer (node) in the sequence.
 */
void BGC(int x)
{
    nbAppel++; // Track the number of recursive calls (search effort)

    // --- CASE 1: Full Sequence Found ---
    if (s >= n2)
    {
        if (s > smax)
            smax++;
        listeCodeBGC.push_back(bgc);

        // Log timing and save full sequence to file
        PE.stop();
        std::ofstream tmp("time.txt", std::ios::app);
        tmp << "\n"
            << PE.seconds();
        tmp.close();

        std::ofstream otp("codes.txt", std::ios::app);
        otp << "\n";
        for (int i : bgc)
        {
            otp << i;
            otp << " ";
        }
        otp.close();
        return;
    }

    // --- CASE 2: Partial Progress Tracking ---
    if (s >= smax)
    {
        // Checkpointing: save the best partial sequence found so far
        if (s > smax)
            smax++;

        PE.stop();

        // ... (File I/O logic for smax progress)
        std::ofstream tmp("time.txt", std::ios::app);
        tmp << "\n"
            << smax << "\n"
            << PE.seconds();
        tmp.close();
        // On stocke dans un fichier texte
        std::ofstream otp("codes.txt", std::ios::app);
        otp << "\n"
            << smax << "\n";
        for (int i : bgc)
        {
            otp << i;
            otp << " ";
        }
        otp.close();
    }

    // --- CASE 3: Pruning by Beckett Constraint ---
    // If a bit has been at '1' for longer than nbMax[n], this branch is invalid.
    if (s > 1 && s - oldS[0] > nbMax[n])
        return;

    // --- CASE 4: Connectivity Analysis (Isolated Vertices) ---
    int cmpt = 0; // Counts neighbors with degree 1
    int xn = 0;   // Candidate for next node
    int bit = 0;  // Index of the flipped bit

    for (int i = 0; i < n; i++)
    {
        int fxi = flip[x][i];
        if (nombreVoisins[fxi] == 1 && avail[fxi])
        {
            cmpt++;
            // Pruning: if more than one neighbor has degree 1, a Hamiltonian path is impossible
            if (cmpt > 1)
                return;

            xn = fxi;
            bit = n - i - 1;
        }
    }

    // --- CASE 5: Mandatory Move (Degree-1 Neighbor) ---
    if (cmpt == 1 && (bi[xn][bit] == 1 || bit == old[0]))
    {
        // Backup state for backtracking
        int oldS0 = oldS[0];

        // Apply Beckett FIFO logic:
        if (bi[xn][bit] == 1)
        { // 0 -> 1: Push bit to the end of the queue
            old.push_back(bit);
            oldS.push_back(s);
        }
        else
        { // 1 -> 0: Only the oldest '1' bit can be flipped back (FIFO)
            old.erase(old.begin());
            oldS.erase(oldS.begin());
        }

        // Move to neighbor xn
        bgc.push_back(xn);
        for (int j = 0; j < n; j++)
            nombreVoisins[flip[xn][j]]--;
        avail[xn] = false;
        s++;

        // Update layer-to-layer edge count (listeVoisinV)
        if (V[x] > V[xn])
            listeVoisinV[V[xn]]--;
        else
            listeVoisinV[V[x]]--;

        BGC(xn); // RECURSION

        // --- BACKTRACK: Revert all changes ---
        s--;
        if (V[x] > V[xn])
            listeVoisinV[V[xn]]++;
        else
            listeVoisinV[V[x]]++;

        for (int j = 0; j < n; j++)
            nombreVoisins[flip[xn][j]]++;

        if (bi[xn][bit] == 1)
        {
            old.pop_back();
            oldS.pop_back();
        }
        else
        {
            old.insert(old.begin(), bit);
            oldS.insert(oldS.begin(), oldS0);
        }
        avail[xn] = true;
        bgc.pop_back();
        return;
    }

    // --- CASE 6: Standard Exploration (Branching) ---
    if (cmpt == 0)
    {
        // Heuristic: Check if there's a "bottleneck" between Hamming layers
        if (!(V[x] == 0 or V[x] == n) && listeVoisinV[V[x] - 1] == 1 && listeVoisinV[V[x]] > 0)
        {
            // Logic similar to mandatory move but specifically for layer bridges...
            // (Repeat logic for bit toggling and backtracking)
            for (int i = 0; i < n; i++)
            {
                bit = n - i - 1;
                if (bi[x][bit] == 0)
                {
                    xn = flip[x][i];
                    if (avail[xn] && (bi[xn][bit] == 1 || bit == old[0]))
                    {
                        int oldS0 = oldS[0];
                        if (bi[xn][bit] == 1)
                        {
                            old.push_back(bit);
                            oldS.push_back(s);
                        }
                        else
                        {
                            old.erase(old.begin());
                            oldS.erase(oldS.begin());
                        }

                        avail[xn] = false;
                        bgc.push_back(xn);

                        for (int j = 0; j < n; j++)
                            nombreVoisins[flip[xn][j]]--;

                        s++;

                        if (V[x] > V[xn])
                            listeVoisinV[V[xn]]--;
                        else
                            listeVoisinV[V[x]]--;

                        BGC(xn);

                        if (V[x] > V[xn])
                            listeVoisinV[V[xn]]++;
                        else
                            listeVoisinV[V[x]]++;

                        s--;

                        for (int j = 0; j < n; j++)
                            nombreVoisins[flip[xn][j]]++;

                        if (bi[xn][bit] == 1)
                        {
                            old.pop_back();
                            oldS.pop_back();
                        }
                        else
                        {
                            old.insert(old.begin(), bit);
                            oldS.insert(oldS.begin(), oldS0);
                        }
                        avail[xn] = true;
                        bgc.pop_back();
                    }
                }
            }
            return;
        }
        else
        {
            // General exploration of all available neighbors
            // (Repeat logic for all valid neighbors under Beckett constraints)
            for (int i = 0; i < n; i++)
            {
                xn = flip[x][i];
                bit = n - i - 1;
                if (avail[xn] && (bi[xn][bit] == 1 || bit == old[0]))
                {
                    int oldS0 = (oldS.empty() ? 0 : oldS[0]);

                    if (bi[xn][bit] == 1)
                    {
                        old.push_back(bit);
                        oldS.push_back(s);
                    }
                    else
                    {
                        old.erase(old.begin());
                        oldS.erase(oldS.begin());
                    }
                    avail[xn] = false;
                    bgc.push_back(xn);

                    for (int j = 0; j < n; j++)
                        nombreVoisins[flip[xn][j]]--;

                    s++;

                    if (V[x] > V[xn])
                        listeVoisinV[V[xn]]--;
                    else
                        listeVoisinV[V[x]]--;

                    BGC(xn);

                    if (V[x] > V[xn])
                        listeVoisinV[V[xn]]++;
                    else
                        listeVoisinV[V[x]]++;

                    s--;
                    for (int j = 0; j < n; j++)
                        nombreVoisins[flip[xn][j]]++;

                    if (bi[xn][bit] == 1)
                    {
                        old.pop_back();
                        oldS.pop_back();
                    }
                    else
                    {
                        old.insert(old.begin(), bit);
                        oldS.insert(oldS.begin(), oldS0);
                    }
                    avail[xn] = true;
                    bgc.pop_back();
                }
            }
            return;
        }
        return;
    }
    return;
}

/**
 * Main execution block: Initializes the environment and launches the BGC search.
 */
int main()
{
    // --- Configuration & Constraints ---
    // Max consecutive 1s allowed for each dimension n (Heuristic balance constraints)
    nbMax = {0, 1, 2, 8, 16, 7, 10, 15, 20};
    n = 5;          // Current dimension (number of bits)
    nbAppel = 0;    // Reset recursive call counter
    n2 = pow(2, n); // Total number of nodes in the hypercube (2^n)
    smax = 0;       // Reset progress tracker

    // --- Data Structure Initialization ---
    bgc = {0};         // Start the sequence with the origin node (0)
    listeCodeBGC = {}; // Clear the results container

    // Initialize visitation map (0 is visited, others are available)
    avail.push_back(false);
    for (int i = 0; i < n2 - 1; i++)
        avail.push_back(true);

    // Initialize the degree of each vertex (each node has n neighbors in an n-cube)
    for (int i = 0; i < n2; i++)
        nombreVoisins.push_back(n);

    // --- Precomputing Lookup Tables ---
    bi = listeBinaire(n); // Binary representation table
    flip = listeFlip(n);  // Adjacency/Transition table
    s = 1;                // Current sequence length (starting with node 0)

    for (int i = 0; i < n + 1; i++)
        nbEltEns.push_back(0); // Reset counts for Hamming weights

    old = {};      // Beckett FIFO queue (active bits)
    oldS = {};     // Entry timestamps for active bits
    bmax = 0;      // Reset worst balance record
    bmin = 100000; // Reset best balance record

    V = listeDes1(n);         // Precompute Hamming weight for all integers
    listeVoisinV = pascal(n); // Precompute layer-to-layer connectivity matrix

    // --- Execution & Performance Measurement ---
    PE.start(); // Start hardware/software timers
    BGC(0);     // Launch the recursive backtracking from node 0
    PE.stop();  // Stop timers after search completion

    // --- Result Output ---
    std::cout << " nbr cycles : " << PE.cycles() << std::endl;
    std::cout << " nbr secondes : " << PE.seconds() << std::endl;
    std::cout << " nbr millisecondes : " << PE.milliseconds() << std::endl;
    cout << "Total codes: " << listeCodeBGC.size() << endl;
    cout << "Nombres d'appels: " << nbAppel << endl;

    PE.clear(); // Final reset of performance counters
    return 0;
}