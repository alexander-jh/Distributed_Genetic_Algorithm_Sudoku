// OptimizationProblem.cpp
//
// Source file for algorithm implementation
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include "OptimizationProblem.hpp"

#include <random>

//------------------- OptimizationProblem Member Functions -------------------------//

void OptimizationProblem::generate_state() {
    uint16_t dim = this->initial->state->dim;
    uint16_t cell = this->initial->state->cell;
    // Transform 4 cell x cell matrices into a dim x dim matrix
    for(uint16_t grid = 0; grid < dim; ++grid) {
        for(uint16_t i = 0; i < cell; ++i) {
            for(uint16_t j = 0; j < cell; ++j) {
                // f:[0,3]x[0,1]x[0,1]->[0,3],[0,3]
                // f(g,i_0,j_0) = (row,col)
                uint16_t row = cell*(grid / cell) + i;
                uint16_t col = cell*(grid % cell) + j;
                // Check if element is missing
                uint16_t num = this->initial->state->state[grid].second[i][j];
                // If missing insert random number in range
                if(num == '*') num = range(generator);
                current->board[row][col] = num;
            }
        }
    }
    current->value = evaluate(current);
    current->next = priority_queue<pair<uint16_t, Node *>>();
}

OptimizationProblem::Node::Node(uint16_t dim) {
    this->board = vector<vector<uint16_t>>(dim, vector<uint16_t>(dim, 0));
    this->next = priority_queue<pair<uint16_t, Node *>>();
    this->value = 0;
}

OptimizationProblem::OptimizationProblem(const char *file) {
    // Get initial state object
    this->initial = new SudokuBoard(file);
    // Calculate goal state
    uint16_t dim = initial->state->dim;
    this->current = new Node(dim);
    this->goal = 3*dim*dim;
    this->range = uniform_int_distribution<uint16_t>(1,dim);
    // Generate current state
    generate_state();
}

uint16_t OptimizationProblem::evaluate(Node *state) const {
    uint16_t value = 0;
    uint16_t grid = initial->state->cell;
    uint16_t dim = initial->state->dim;
    vector<unordered_set<uint16_t>> cell = vector<unordered_set<uint16_t>>(dim);
    for(int i = 0; i < dim; i++) {
        unordered_set<uint16_t> row;
        unordered_set<uint16_t> col;
        for(uint16_t j = 0; j < dim; j++) {
            // Count unique in each possible row
            row.insert(state->board[i][j]);
            // Count unique in each possible col
            col.insert(state->board[j][i]);
            // Count unique in each grid g(i,j) for each element in sub-matrix,
            // Inverse of transformation f described in generate_state()
            // f':[0,3],[0,3]->[0,3]x[0,1]x[0,1]
            // f'(i,j) = (g,row,col)
            uint16_t g = grid*(i / grid) + (j / grid);
            cell[g].insert(state->board[i][j]);
        }
        value += row.size() + col.size();
    }
    for(const auto& it : cell) {
        value += it.size();
    }
    return value;
}

void OptimizationProblem::report_results(Node *node, uint16_t goal, uint64_t total) {
    uint16_t dim = node->board.size();
    if(node->value == goal) {
        printf("\nSolution found!\n");
    } else {
        printf("\nNon goal state found.\n");
        printf("\nBEST SOLUTION:\n");
    }
    for(uint16_t i = 0; i < dim; ++i) {
        for(uint16_t j = 0; j < dim; ++j) {
            printf("%d\t", node->board[i][j]);
        }
        printf("\n");
    }
    printf("Iterations:\t%lu\n", total);
}

OptimizationProblem::OptimizationProblem() {
    current = nullptr;
    initial = nullptr;
    goal = 0;
}

//------------------- HillClimber Member Functions -------------------------//

HillClimber::HillClimber(const char *file) {
    state = new OptimizationProblem(file);
}

OptimizationProblem::Node *HillClimber::hill_climb() {
    bool goal_found = false;
    // Counter for no improvement in next candidate state
    uint16_t no_improvement = 0;
    while(!goal_found && state->total_iterations < UINT64_MAX) {
        // Get successors of current node
        goal_found = generate_successors();
        // Report current node and next best candidate
        printf("Current Value: %d\t\t Successor Value: %d\n", state->current->value,
               state->current->next.top().first);
        // If candidate state is better set current state to candidate
        if(state->current->value < state->current->next.top().first || goal_found) {
            state->current = state->current->next.top().second;
        } else {
        // Otherwise increment plateau counter
            no_improvement++;
        // If at plateau reset algorithm
            if(no_improvement == PLATEAU) {
                state->generate_state();
                no_improvement = 0;
            }
        }
        state->total_iterations++;
    }
    return state->current;
}


bool HillClimber::generate_successors() {
    // Adjust a singular value of missing places from current state
    for(auto &it: state->initial->state->missing) {
        uint16_t row = it.first;
        uint16_t col = it.second;
        // Test all possible changes of single position value not equal
        // to the current value
        for (uint16_t i = 1; i < state->initial->state->dim+1; ++i) {
            if(state->current->board[row][col] != i) {
                auto *test = new OptimizationProblem::Node(state->current->board.size());
                test->board = state->current->board;
                uint16_t num = state->range(state->generator);
                // Update value on board
                test->board[row][col] = num;
                test->value = state->evaluate(test);
                state->current->next.push({test->value, test});
                if(test->value == state->goal) return true;
            }
        }
    }
    return false;
}

//------------------- Genetic Member Functions -------------------------//
#ifndef HILL

Genetic::Parent::Parent() {
    this->value = 0;
}

Genetic::Genetic(const char *file) {
    // Initialize mutexes
    pthread_mutex_init(&mutex_pop_in, nullptr);
    pthread_mutex_init(&mutex_pop_out, nullptr);
    // Initialize attr
    pthread_attr_init(&attr);
    auto *start = new OptimizationProblem(file);
    // Create initial test state
    for(const auto& it : start->current->board)
        initial.push_back(it);
    range = uniform_int_distribution<uint16_t>(1,start->initial->state->dim);
    // Move all missing pair locations over
    for(auto it : start->initial->state->missing) {
        missing.push_back(it);
    }
    goal = start->goal;
    delete start;
    iterations = 0;
    reset_search(0);
}

void Genetic::reset_search(uint16_t elite) {
    for(int i = 0; i < POP_MODIFIER*MAX_THREADS - elite; ++i) {
        auto *p = new Parent;
        for(int j = 0; j < missing.size(); ++j) p->genes.push_back(range(generator));
        population.push_back(p);
    }
    // Generate initial fitness for each parent
    buffer_in = 0;
    for(auto & i : thread_pool)
        pthread_create(&i, &attr, &Genetic::update_help, this);
    for(auto & i : thread_pool)
        pthread_join(i, nullptr);
}

bool Genetic::genetic_run() {
    uniform_int_distribution<uint16_t> rand;
    rand = uniform_int_distribution<uint16_t>(1,10);
    goal_met = false;
    // Recombine until goal is found
    while(!goal_met) {
        if(rand(generator) == 1) {
            printf("\n RESET SEARCH \n");
            vector<Parent *> elite(POP_MODIFIER);
            sort(population.begin(), population.end(), Parent::comparator());
            for(int i = 0; i < POP_MODIFIER; ++i)
                    elite[i] = new Parent(population.back());
            population.clear();
            reset_search(POP_MODIFIER);
            while(!elite.empty()) {
                population.push_back(new Parent(elite.back()));
                elite.pop_back();
            }
        }
        for(auto & i : thread_pool)
            pthread_create(&i, &attr, &Genetic::thread_help, this);
        for(auto & i : thread_pool)
            pthread_join(i, nullptr);
        iterations++;
    }
    return goal_met;
}

void *Genetic::thread_run() {
    // Obtain mutex lock to take the two largest parents
    vector<Parent *> a;
    pthread_mutex_lock(&this->mutex_pop_out);
    // Perform selection
    if(!this->goal_met) selection(a);
    pthread_mutex_unlock(&this->mutex_pop_out);
    if(this->goal_met) pthread_exit(this);
    // Crossover
    a = crossover(a[0],a[1]);
    printf("Parent A Value: %d\tParent B Value: %d\n",a[0]->value, a[1]->value);
    // Reacquire lock
    pthread_mutex_lock(&this->mutex_pop_out);
    if(!this->goal_met) {
        this->goal_met = (a[0]->value == this->goal || a[1]->value == this->goal);
        this->population.push_back(a[0]);
        this->population.push_back(a[1]);
    }
    pthread_mutex_unlock(&this->mutex_pop_out);
    pthread_exit(this);
}

void Genetic::genetic_report(bool status) {
    int count = 0;
    if(status) {
    // Insert best candidate into initial frame
        Parent *max;
        auto pop = population.begin();
        while((*pop)->value != goal) {
            pop++;
        }
        max = *pop;
        for(auto it : missing) {
            initial[it.first][it.second] = max->genes[count];
            count++;
        }
        printf("\nSOLUTION FOUND!\n");
    // Print results
        for(auto & i : initial) {
            for(int j = 0; j < initial.size(); ++j) {
                printf("%d\t", i[j]);
            }
            printf("\n");
        }
        printf("Iterations: %lu\n", iterations);
    } else {
        printf("\nNO SOLUTION FOUND!\n");
    }
}

vector<Genetic::Parent *> Genetic::crossover(Genetic::Parent *a, Genetic::Parent *b) {
    auto *c = new Parent, *d = new Parent;
    uniform_real_distribution<double> prob(0,1);
    for(int i = 0; i < a->genes.size(); ++i) {
        double roll = prob(generator);
        (roll <= 0.5) ? c->genes.push_back(a->genes[i]) : c->genes.push_back(b->genes[i]);
        roll = prob(generator);
        (roll <= 0.5) ? d->genes.push_back(a->genes[i]) : d->genes.push_back(b->genes[i]);
    }
    // One percent chance of mutation
    if(this->mutate(this->generator) == 1) mutation(c);
    if(this->mutate(this->generator) == 1) mutation(d);
    // Get new fitness value
    c->value = fitness(c);
    d->value = fitness(d);
    // Get fist max of the four parents
    vector<pair<uint16_t, Parent *>> candidates = {{a->value,a},{b->value,b},
                                                   {c->value,c},{d->value,d}};
    sort(candidates.begin(),candidates.end());
    vector<Parent *> ret;
    for(uint16_t i = 0; i < 4; ++i) {
        if(i < 2) {
            ret.push_back(candidates.back().second);
        } else {
            delete(candidates.back().second);
        }
        candidates.pop_back();
    }
    return ret;
}

void Genetic::selection(vector<Parent *> &pool) {
    uint16_t loc = uniform_int_distribution<uint16_t>(0, population.size() - 1)(generator);
    pool.push_back(population[loc]);
    population.erase(population.begin() + loc);
    pool.push_back(population.back());
    population.pop_back();
}

void Genetic::mutation(Genetic::Parent *a) {
    uniform_int_distribution<uint16_t> dis;
    dis = uniform_int_distribution<uint16_t>(0,a->genes.size()-1);
    a->genes[dis(generator)] = range(generator);
}

void *Genetic::update_fitness(void *arg) {
    vector<uint16_t> pos;
    // Lock the population vector
    while(buffer_in < this->population.size()) {
        pthread_mutex_lock(&this->mutex_pop_out);
        // Remove last individual
        pos.push_back(buffer_in);
        buffer_in++;
        // Drop lock and calculate fitness
        pthread_mutex_unlock(&this->mutex_pop_out);
        if(buffer_in < this->population.size())
            this->population[pos.back()]->value = fitness(this->population[pos.back()]);
    }
    pthread_exit(this);
}

uint16_t Genetic::fitness(Genetic::Parent *a) {
    uint16_t value = 0;
    uint16_t cell = Board::isqrt(initial.size());
    uint16_t dim = initial.size();
    vector<unordered_set<uint16_t>> grid = vector<unordered_set<uint16_t>>(dim);
    auto *temp = new vector<vector<uint16_t>>(initial);
    // Insert missing elements
    uint16_t index = 0;
    for(const auto& it : missing) {
        (*temp)[it.first][it.second] = a->genes[index];
        index++;
    }
    for(int i = 0; i < dim; i++) {
        unordered_set<uint16_t> row;
        unordered_set<uint16_t> col;
        for(uint16_t j = 0; j < dim; j++) {
            // Count unique in each possible row
            row.insert((*temp)[i][j]);
            // Count unique in each possible col
            col.insert((*temp)[j][i]);
            // Count unique in each grid g(i,j) for each element in sub-matrix,
            // Inverse of transformation f described in generate_state()
            // f':[0,3],[0,3]->[0,3]x[0,1]x[0,1]
            // f'(i,j) = (g,row,col)
            uint16_t g = cell * (i / cell) + (j / cell);
            grid[g].insert((*temp)[i][j]);
        }
        value += row.size() + col.size();
    }
    for(const auto& it : grid) {
        value += it.size();
    }
    delete temp;
    return value;
}
#endif