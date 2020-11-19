// OptimizationProblem.hpp
//
// Header files for optimization tools, hill climbing algorithm, and
// the genetic algorithm.
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include "SudokuBoard.hpp"
#include <pthread.h>

#define MAX_THREADS     16

// Simple XOR hash algorithm for dealing with undefined objects
struct pair_hash {
    template <class T1, class T2>
    size_t operator () (std::pair<T1, T2> const &pair) const {
        size_t h1 = std::hash<T1>()(pair.first);
        size_t h2 = std::hash<T2>()(pair.second);
        return h1 ^ h2;
    }
};

// Generic tool kit for translating the initial SudokuBoard representation
// into a full matrix for easier parsing
class OptimizationProblem {
public:
    // Member struct to represent current and initial public states
    // @param board     -   NxN matrix to represent the game state
    // @param next      -   Priority to represent candidate states for the
    //                      hill climb algorithm
    // @param value     -   Heuristic value for the present state
    typedef struct Node{
        vector<vector<uint16_t>>                   board;
        priority_queue<pair<uint16_t, Node*>>      next;
        uint16_t                                   value;
        // Explicit constructor to allocate memory space
        // @param dim   -   overall dimension of the board matrix
        explicit Node(uint16_t dim);
        // Overloaded copy constructor to verify pointer objects are
        // copied over
        // @param node  -   Node to be duplicated
        inline Node& operator= (const Node& node) {
            board = node.board;
            for(auto &it : node.board)
                board.emplace_back(it);
            next = priority_queue<pair<uint16_t, Node*>>();
            value = node.value;
            return *this;
        }
    }Node;
    // Current state representation
    Node *current;
    // Number of iterations
    uint64_t iterations = 0;
    // Default constructor
    OptimizationProblem();
    // Constructor with file input
    // @param *file -   pointer to file location
    explicit OptimizationProblem(const char *file);
    // Builds initial state from SudokuBoard whenever local max is hit or
    // for initial construction
    void generate_state();
    // Runs evaluation function on current board state. Heuristic evaluates the
    // state by counting the sum total of all unique elements in each row, column,
    // and in each sub-matrix. Higher value implies better state of the board.
    // @param *state    -   current state to be evaluated
    uint16_t evaluate(Node *state);
    // Reports the results from the algorithm
    // @param status    -   boolean to distinguish if solution was found
    void report_results(bool status) const;
    // Initial board representation from parsing class
    SudokuBoard *initial;
    // Solution state constant, goal = 3 * dim
    uint16_t goal;
    // Range of random number such that [1,dim+1]
    uniform_int_distribution<uint16_t> range;
    // Psuedo-random number generator, more reliable than rand()
    random_device generator;
};

// Implementation of the Hill Climbing Algorithm
class HillClimber {
public:
    // Constructor for when file is passed
    // @param *file -   pointer to file location
    explicit HillClimber(const char *file);
    // Constructor for when no file is passed. Default behavior is to
    // terminate the program
    inline HillClimber(){SudokuBoard();};
    // Default destructor
    ~HillClimber() = default;
    // Control flow for the algorithm
    bool hill_climb();
    // Reports the results
    inline void hill_report(bool status) {state->report_results(status);};
private:
    // Present state representation
    OptimizationProblem *state;
    // Fill the scores queue in Node for all possible successors from
    // changing a single value
    bool generate_successors();
};

// Distributed implementation of the Genetic algorithm using p_threads.
//
// Heuristic is the same as in Hill Climbing algorithm. Genome is represented
// as a 1-dimensional array containing the originally empty elements.
//
// Uses a simple random sliding crossover
//
// Selection is tournament style. Selecting 4 random candidates from the
// population, then selecting the best two for breeding
//
// Mutation is a 1% chance, randomly selects and changes the value with
// the highest count in the current state
//
// Resets are done every 100 iterations if solution not found
//
// Threads are being used based on a cell implementation. Mutual exclusion is
// help on the population whenever selection or reinsertion is done. However,
// the total population is 2*THREAD_MAX with each thread only holding 2 parents.
// With a fluctuating population between runs this increases randomness of the
// selection process.
class Genetic {
public:
    // Parent 1-d array representation
    // @param genes -   1-d array of the genome of a given state
    // @param value -   current heuristic value of the state
    typedef struct Parent {
        vector<uint16_t>                           genes;
        uint16_t                                   value;
        // Default no arg constructor
        explicit Parent();
        // Copy constructor for memory protections
        inline explicit Parent(Parent *node) {
            genes = vector<uint16_t>(node->genes.size());
            for(int i = 0; i < node->genes.size(); ++i) {
                genes[i] = node->genes[i];
            }
            value = node->value;
        }
        // Copy constructor for operator =
        inline Parent& operator= (const Parent& node) {
            genes = vector<uint16_t>(node.genes.size());
            for(int i = 0; i < node.genes.size(); ++i) {
                genes[i] = node.genes[i];
            }
            value = node.value;
            return *this;
        }
    } Parent;
    // Constructor for when file is passed
    // @param *file -   pointer to file location
    explicit Genetic(const char *file);
    // Constructor for when no file is passed. Default behavior kills program.
    inline Genetic(){SudokuBoard();};
    // Default destructor
    ~Genetic() = default;
    // Runs the algorithm
    bool genetic_run();
    // Reports the final state
    // @param status    -   boolean signifying if goal state achieved
    void genetic_report(bool status);
private:
    // Iteration and reset count for the program
    uint16_t iterations, resets;
    // Population mutex lock for initial input
    pthread_mutex_t mutex_pop_in;
    // Population mutex lock for taking out individuals
    pthread_mutex_t mutex_pop_out;
    // Default pthread attr
    pthread_attr_t attr;
    // Array for pthread pool
    pthread_t thread_pool[MAX_THREADS];
    // Direct (i,j) position of each missing element in original state
    vector<pair<uint16_t, uint16_t>> missing;
    // Goal state 3*dim
    uint16_t goal;
    // End state achieved declared in the private name space to message
    // threads to quit as soon as a goal state is found
    bool goal_met;
    // Boolean declared at name space for initial construction of the population.
    // Stalls threads until population is empty to prevent race conditions
    bool thread_done;
    // Range of random number such that [1,dim+1]
    uniform_int_distribution<uint16_t> range;
    // Generate random number
    random_device generator;
    // Mutation chance generator
    uniform_int_distribution<uint16_t> mutate = uniform_int_distribution<uint16_t>(1,100);
    // NxN representation of the board state. Used for copy to calculate heuristic
    vector<vector<uint16_t>> initial;
    // Population of parent nodes
    vector<Parent *> population;
    // Hash set for ensuring proper ordering of the gene string
    unordered_set<pair<uint16_t, uint16_t>, pair_hash> miss;
    // Conducts a random sliding window crossover. Randomly generates an index position
    // in [0,gene.size]. Direct copies from first parent, then fills rest from second
    // parent. Second child is created doing the opposite.
    // @param *a    -   first candidate parent node
    // @param *b    -   second candidate parent node
    void crossover(Parent *a, Parent *b);
    // Conducts mutation, 1 in 100 chance of mutation
    // @param *a    -   candidate node for mutation
    void mutation(Parent *a);
    // Evaluation function, the heuristic algorithm parallels the hill climbing implementation.
    // Had to be recreated to load the current gene array into the board prior to evaluation
    // @param *a    -   node to be evaluated
    uint16_t fitness(Parent *a);
    // Selection is conducted tournament style. 4 candidates are chosen randomly from population.
    // The two "fittest" candidates are removed from the population and returned to the calling
    // thread
    void selection(vector<Parent *> &pool);
    // Called from the first pthread creation where candidate states are evaluated. Necessary to
    // ensure program is free from race conditions. Simply removes parents, calculates fitness
    // and returns them to the population. Is distributed to increase efficiency.
    void *update_fitness(void *arg);
    // Helper method for update_fitness. Pthreads require a static function declaration so
    // the helper is needed.
    inline static void *update_help(void *context) {return ((Genetic *)context)->update_fitness(nullptr);};
    // Thread function for conducting selection, conducting crossover, mutation, evaluation,
    // and reinsertion to population.
    //
    // Selection is in the critical zone and is protected by mutex locks. Crossover and evaluation
    // are done in parallel (or concurrently). Reinsertion to population is done with mutex locks.
    void *thread_run();
    // Similar helper method as previously detailed.
    inline static void *thread_help(void *context) {return ((Genetic *) context)->thread_run();};
    // Constructs the initial game state, or resets the state if stuck in a min or max.
    void reset_search();
};
