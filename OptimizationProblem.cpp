// OptimizationProblem.cpp
//
// Source file for algorithm implementation
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include "OptimizationProblem.hpp"

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

uint16_t OptimizationProblem::evaluate(Node *state) {
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

void OptimizationProblem::report_results(bool status) const {
    if(status) {
        printf("\nSolution found!\n");
        for(uint16_t i = 0; i < initial->state->dim; ++i) {
            for(uint16_t j = 0; j < initial->state->dim; ++j) {
                printf("%d\t", current->board[i][j]);
            }
            printf("\n");
        }
        printf("Iterations:\t%lu\n",iterations);
    } else {
        printf("\nNo solution found.\n");
    }
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

bool HillClimber::hill_climb() {
    bool goal_found = false;
    OptimizationProblem::Node *curr;
    while(!goal_found) {
        curr = state->current;
        goal_found = generate_successors();
        printf("Current Value: %d\t\t Successor Value: %d\n", curr->value, curr->next.top().first);
        if(curr->value < curr->next.top().first || goal_found) {
            state->current = curr->next.top().second;
        } else {
            printf("No valid successor reset game state.\n");
            state->generate_state();
        }
        if (state->iterations != UINT64_MAX) {state->iterations++;} else {break;}
    }
    return goal_found;
}


bool HillClimber::generate_successors() {
    // Adjust value of missing places
    for(auto &it: state->initial->state->missing) {
        auto *test = new OptimizationProblem::Node(state->current->board.size());
        uint16_t row = it->first;
        uint16_t col = it->second;
        test->board = state->current->board;
        uint16_t num = state->range(state->generator);
        // Update value on board
        test->board[row][col] = num;
        test->value = state->evaluate(test);
        state->current->next.push({test->value, test});
        if(state->current->next.top().first == state->goal) return true;
    }
    return false;
}

//------------------- Genetic Member Functions -------------------------//

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
    for(const auto it : start->current->board)
        initial.push_back(it);
    range = uniform_int_distribution<uint16_t>(1,start->initial->state->dim);
    // Move all missing pair locations over
    for(auto it : start->initial->state->missing) {
        miss.insert(*it);
    }
    goal = start->goal;
    delete start;
    resets = 0;
    reset_search();
}

void Genetic::reset_search() {
    // Create initial population with 2*MAX_THREADS generated children
    missing.clear();
    uint16_t dim = initial.size();
    uint16_t cell = Board::isqrt(dim);
    for(uint16_t i = 0; i < 4*MAX_THREADS; ++i) {
        // Generate children by mapping random gene not already in sub-matrix at
        // missing positions
        auto *p = new Parent;
        for(uint16_t g = 0; g < dim; ++g) {
            unordered_set<uint16_t> cell_count;
            // Get existing elements from each sub-matrix
            for(uint16_t j = 0; j < cell; ++j) {
                for (uint16_t k = 0; k < cell; ++k) {
                    uint16_t R = cell * (g / cell) + j;
                    uint16_t C = cell * (g % cell) + k;
                    auto it = miss.find(make_pair(R,C));
                    if(it == miss.end()) {
                        cell_count.insert(initial[R][C]);
                    } else if (i == 0) {
                        missing.push_back(*it);
                    }
                }
            }
            // Insert unique missing genes
            while(cell_count.size() != dim) {
                uint16_t val = range(generator);
                while(cell_count.find(val) != cell_count.end()) {
                    val = range(generator);
                }
                p->genes.push_back(val);
                cell_count.insert(val);
            }
        }
        population.push_back(p);
    }
    // Generate initial fitness for each parent
    for(int i = 0; i < MAX_THREADS; ++i)
        pthread_create(&thread_pool[i], &attr, &Genetic::update_help, this);
    // Spinlock until priority queue has emptied
    while(!(thread_done = population.empty()));
    for(int i = 0; i < MAX_THREADS; ++i)
        pthread_join(thread_pool[i], nullptr);
    iterations = 0;
}

bool Genetic::genetic_run() {
    goal_met = false;
    // Recombine until goal is found
    while(!goal_met) {
        if(iterations == 100) {
            resets++;
            printf("\n RESET SEARCH \n");
            population.clear();
            reset_search();
        }
        thread_done = false;
        for(int i = 0; i < MAX_THREADS; ++i)
            pthread_create(&thread_pool[i], &attr, &Genetic::thread_help, this);
        for(int i = 0; i < MAX_THREADS; ++i)
            pthread_join(thread_pool[i], nullptr);
        iterations++;
    }
    return goal_met;
}

void *Genetic::thread_run() {
    // Obtain mutex lock to take the two largest parents
    vector<Parent *> a;
    pthread_mutex_lock(&this->mutex_pop_out);
    if(!this->goal_met) {
        // Perform selection
        selection(a);
    }
    pthread_mutex_unlock(&this->mutex_pop_out);
    if(this->goal_met) pthread_exit(nullptr);
    // Crossover
    crossover(a[0],a[1]);
    // One percent chance of mutation
    if(this->mutate(this->generator) == 1) mutation(a[0]);
    if(this->mutate(this->generator) == 1) mutation(a[1]);
    // Get new fitness value
    a[0]->value = fitness(a[0]);
    a[1]->value = fitness(a[1]);
    printf("Parent A Value: %d\tParent B Value: %d\n",a[0]->value, a[1]->value);
    // Reacquire lock
    pthread_mutex_lock(&this->mutex_pop_out);
    if(!this->goal_met) {
        this->goal_met = (a[0]->value == this->goal || a[1]->value == this->goal);
        this->population.push_back(a[0]);
        this->population.push_back(a[1]);
    }
    pthread_mutex_unlock(&this->mutex_pop_out);
    pthread_exit(nullptr);
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
        for(int i = 0; i < initial.size(); ++i) {
            for(int j = 0; j < initial.size(); ++j) {
                printf("%d\t", initial[i][j]);
            }
            printf("\n");
        }
        printf("Iterations: %d\n", iterations+(100*resets));
    } else {
        printf("\nNO SOLUTION FOUND!\n");
    }
}

void Genetic::crossover(Genetic::Parent *a, Genetic::Parent *b) {
    uint16_t mid = range(generator);
    auto *c = new Parent, *d = new Parent;
    for(uint16_t i = 0; i < mid; i++) {
        c->genes.push_back(b->genes[i]);
        d->genes.push_back(a->genes[i]);
    }
    for(int i = mid; i < a->genes.size(); ++i) {
        c->genes.push_back(a->genes[i]);
        d->genes.push_back(b->genes[i]);
    }
    a = new Parent(c), b = new Parent(d);
    delete c, d;
}

// Tournament based selection
void Genetic::selection(vector<Parent *> &pool) {
    // Randomly select 4 parents
    for(uint16_t i = 0; i < 4; ++i) {
        uniform_int_distribution<uint16_t> pos =
                uniform_int_distribution<uint16_t> (0, population.size()-1);
        uint16_t loc = pos(generator);
        if(i < 2) {
            pool.push_back(new Parent(population[loc]));
            population.erase(population.begin()+loc);
        } else {
            uint16_t min;
            (pool[0]->value >= pool[1]->value) ? (min = 1) : (min = 0);
            if(population[loc]->value > pool[min]->value) {
                population.push_back(new Parent(pool[min]));
                pool[min] = new Parent(population[loc]);
                population.erase(population.begin()+loc);
            }
        }
    }
}

void Genetic::mutation(Genetic::Parent *a) {
    unordered_map<uint16_t, uint16_t> count;
    uint16_t max, ele;
    for(auto it : a->genes) {
        auto pos = count.find(it);
        if(pos == count.end()) {
            count.insert({it, 1});
            if(max == 0) {
                max = 1;
                ele = it;
            }
        } else {
            pos->second++;
            if(pos->second > max) {
                max = pos->second;
                ele = pos->first;
            }
        }
    }
    uniform_int_distribution<uint16_t> dis = uniform_int_distribution<uint16_t>(1,max);
    uint16_t place = dis(generator);
    for(auto it : a->genes) {
        if(it == ele) {
            place--;
            if(place == 0) it = range(generator);
        }
    }
}

void *Genetic::update_fitness(void *arg) {
    // Since 2*MAX_THREAD is population size each thread takes two
    Genetic::Parent *a, *b, *c, *d;
    // Lock the population vector
    pthread_mutex_lock(&this->mutex_pop_out);
    // Remove last individual
    a = this->population.back();
    this->population.pop_back();
    b = this->population.back();
    this->population.pop_back();
    c = this->population.back();
    this->population.pop_back();
    d = this->population.back();
    this->population.pop_back();
    pthread_mutex_unlock(&this->mutex_pop_out);
    // Drop lock and calculate fitness
    a->value = fitness(a);
    b->value = fitness(b);
    c->value = fitness(c);
    d->value = fitness(d);
    // Re-obtain lock to insert updated individual
    pthread_mutex_lock(&this->mutex_pop_in);
    while(!this->thread_done);
    this->population.push_back(a);
    this->population.push_back(b);
    this->population.push_back(c);
    this->population.push_back(d);
    pthread_mutex_unlock(&this->mutex_pop_in);
    pthread_exit(this);
}

uint16_t Genetic::fitness(Genetic::Parent *a) {
    uint16_t value = 0;
    uint16_t cell = Board::isqrt(initial.size());
    uint16_t dim = initial.size();
    vector<unordered_set<uint16_t>> grid = vector<unordered_set<uint16_t>>(dim);
    vector<vector<uint16_t>> temp = *new vector<vector<uint16_t>>(initial);
    // Insert missing elements
    uint16_t index = 0;
    for(const auto& it : missing) {
        temp[it.first][it.second] = a->genes[index];
        index++;
    }
    for(int i = 0; i < dim; i++) {
        unordered_set<uint16_t> row;
        unordered_set<uint16_t> col;
        for(uint16_t j = 0; j < dim; j++) {
            // Count unique in each possible row
            row.insert(temp[i][j]);
            // Count unique in each possible col
            col.insert(temp[j][i]);
            // Count unique in each grid g(i,j) for each element in sub-matrix,
            // Inverse of transformation f described in generate_state()
            // f':[0,3],[0,3]->[0,3]x[0,1]x[0,1]
            // f'(i,j) = (g,row,col)
            uint16_t g = cell * (i / cell) + (j / cell);
            grid[g].insert(temp[i][j]);
        }
        value += row.size() + col.size();
    }
    for(const auto& it : grid) {
        value += it.size();
    }
    return value;
}
