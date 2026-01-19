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
vector<int> p2;            // Vecteur des puissances de 2 plus petites que n**2
vector<bool> p2Bool;       // Vecteur à true pour des puissances de 2 plus petites que n**2, false sinon
vector<vector<int>> bi;    // Precomputed binary representation for each integer < 2^n
vector<vector<int>> flip;  // Adjacency matrix: flip[i][j] is the neighbor of i by toggling bit j
vector<int> nombreVoisins; // Current degree of each vertex (number of unvisited neighbors)
vector<int> V;             // Hamming weight: V[i] stores the number of set bits (1s) in i
vector<int> nbEltEns;      // Population count: number of remaining elements for each Hamming weight

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
 * Generates a list of all powers of two within the 2^n range.
 * * These values correspond to integers with exactly one bit set (Hamming weight = 1).
 * * @param n The number of bits (dimension).
 * @return A vector containing all powers of two from 2^0 to 2^(n-1).
 */
vector<int> listepower2(int n)
{
    vector<int> res;
    for (int i = 0; i < pow(2, n); i++)
    {
        // Use the previously defined bitwise check
        if (pw2(i))
        {
            res.push_back(i);
        }
    }
    return res;
}

/**
 * Creates a boolean lookup table (mask) to identify powers of two.
 * Each index i in the returned vector is True if i is a power of two,
 * and False otherwise.
 * @param n The number of bits defining the range [0, 2^n - 1].
 * @return A vector of booleans of size 2^n.
 */
vector<bool> listepower2Bool(int n)
{
    vector<bool> res;
    for (int i = 0; i < pow(2, n); i++)
    {
        // Check if the current index is a power of two using bitwise logic
        if (pw2(i))
        {
            res.push_back(true);
        }
        else
        {
            res.push_back(false);
        }
    }
    return res;
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
                cmpt++;
            else
            {
                // Sequence of 1s broken, update the max if necessary
                if (cmpt > max)
                    max = cmpt;

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
 * Recursive Backtracking function for Beckett Gray Code discovery.
 * Features a specific pruning check for Power-of-Two nodes.
 */
void BGC(int x)
{
    nbAppel++;

    // --- 1. SUCCESS CONDITION ---
    if (s >= n2)
    {
        if (s > smax)
            smax++;

        // Log timing and save full sequence
        listeCodeBGC.push_back(bgc);

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

    // --- 2. PROGRESS LOGGING ---
    if (s >= smax)
    {
        if (s > smax)
            smax++;

        // Periodic save of the longest partial path found so far
        PE.stop();
        std::ofstream tmp("time.txt", std::ios::app);
        tmp << "\n"
            << smax << "\n"
            << PE.seconds();
        tmp.close();

        // ... (File I/O for partial results)

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

    // --- 3. BECKETT CONSTRAINT CHECK ---
    // FIFO rule: If a bit stays at '1' too long, the code is no longer "Beckett"
    if (s > 1 && s - oldS[0] > nbMax[n])
        return;

    // --- 4. POWER-OF-TWO PRUNING (Unique to this version) ---
    // If we are currently on a Power of Two, verify if other Power-of-Two
    // nodes are still reachable. If not, this branch is a dead end.
    if (p2Bool[x])
    {
        bool flag = false;
        for (int i : p2)
        {
            if (avail[i])
            {
                flag = true;
                break;
            }
        }
        if (!flag)
            return;
    }

    // --- 5. TOPOLOGY ANALYSIS: ISOLATED VERTICES ---
    int cmpt = 0; // Number of neighbors with degree 1
    int xn = 0;
    int bit = 0;

    for (int i = 0; i < n; i++)
    {
        int fxi = flip[x][i];
        if (nombreVoisins[fxi] == 1 && avail[fxi])
        {
            cmpt++;

            // Dead-end detection: A Hamiltonian path can only handle one degree-1 node

            if (cmpt > 1)
                return;

            xn = fxi;
            bit = n - i - 1;
        }
    }

    // --- 6. FORCED MOVE (One isolated neighbor exists) ---
    if (cmpt == 1 && (bi[xn][bit] == 1 || bit == old[0]))
    {
        int oldS0 = oldS[0];

        // Manage Beckett FIFO Queue
        if (bi[xn][bit] == 1)
        { // Adding a '1'
            old.push_back(bit);
            oldS.push_back(s);
        }
        else
        { // Removing a '1' (Must be the oldest one)
            old.erase(old.begin());
            oldS.erase(oldS.begin());
        }

        // Update state and recurse
        bgc.push_back(xn);

        for (int j = 0; j < n; j++)
            nombreVoisins[flip[xn][j]]--;

        avail[xn] = false;
        s++;
        BGC(xn);

        // --- BACKTRACK: Restore state ---
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
        return;
    }

    // --- 7. BRANCHING EXPLORATION (No isolated neighbors) ---
    if (cmpt == 0)
    {
        for (int i = 0; i < n; i++)
        {
            xn = flip[x][i];
            bit = n - i - 1;

            // Check availability and Beckett FIFO validity
            if (avail[xn] && (bi[xn][bit] == 1 || bit == old[0]))
            {
                int oldS0 = (oldS.empty() ? 0 : oldS[0]);

                // Similar state management to Case 6...
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

                bgc.push_back(xn);

                for (int j = 0; j < n; j++)
                    nombreVoisins[flip[xn][j]]--;

                s++;
                avail[xn] = false;

                BGC(xn); // Recursive call

                // Backtrack
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

/**
 * Main execution block: Sets up the Beckett Gray Code search environment.
 * This version focuses on Power-of-Two tracking for pruning.
 */
int main()
{
    // --- 1. CONFIGURATION ---
    // Maximum allowed run-length of 1s for each bit position
    nbMax = {0, 1, 2, 8, 16, 7, 10, 15, 20};
    n = 5;          // Dimension of the hypercube
    nbAppel = 0;    // Counter for recursive calls
    n2 = pow(2, n); // Total nodes (32 for n=5)
    smax = 0;       // High-water mark for sequence length found

    // --- 2. SEQUENCE INITIALIZATION ---
    bgc = {0};         // Sequence starts at the origin (00000)
    listeCodeBGC = {}; // Reset results container

    // --- 3. GRAPH STATE SETUP ---
    // visitation map: 0 is visited (false), all other 2^n - 1 nodes are available (true)
    avail.push_back(false);
    for (int i = 0; i < n2 - 1; i++)
        avail.push_back(true);

    // Initial degree of each vertex in an n-dimensional hypercube
    for (int i = 0; i < n2; i++)
        nombreVoisins.push_back(n);

    // --- 4. PRECOMPUTATION (Lookup Tables) ---
    bi = listeBinaire(n);        // Binary patterns table
    flip = listeFlip(n);         // Adjacency/Neighbor table
    p2 = listepower2(n);         // List of integers that are powers of 2
    p2Bool = listepower2Bool(n); // Boolean mask for O(1) power-of-two checks

    // --- 5. SEARCH STATE RESET ---
    s = 1;     // Current path length (only '0' is in)
    old = {};  // Beckett FIFO queue (bit indices)
    oldS = {}; // Timestamps for bit entry
    bmax = 0;
    bmin = 100000;

    // --- 6. EXECUTION & BENCHMARKING ---
    PE.start(); // Start high-resolution hardware timers
    BGC(0);     // Start recursive backtracking search
    PE.stop();  // Stop timers

    // --- 7. OUTPUT STATISTICS ---
    std::cout << " nbr cycles : " << PE.cycles() << std::endl;
    std::cout << " nbr secondes : " << PE.seconds() << std::endl;
    std::cout << " nbr millisecondes : " << PE.milliseconds() << std::endl;
    cout << "Total codes: " << listeCodeBGC.size() << endl;
    cout << "Nombres d'appels: " << nbAppel << endl;

    PE.clear(); // Reset performance class
    return 0;
}