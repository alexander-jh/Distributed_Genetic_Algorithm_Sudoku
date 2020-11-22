// SudokuBoard.cpp
//
// Source implementation for the game board representation.
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include "SudokuBoard.hpp"

Board::Board(uint16_t dim) {
    this->dim = dim;
    this->cell = isqrt(dim);
    // Allocate heap space for missing
    this->missing = vector<pair<uint16_t, uint16_t>>();
    // Create individual pre-sized vectors
    for(uint16_t i = 0; i < dim; i++)
    this->state.emplace_back(*new unordered_set<uint16_t>,
                             vector<vector<uint16_t>>(this->cell, vector<uint16_t>(this->cell, 0)));
}

Board::~Board() {
    delete &missing;
    delete &state;
}

SudokuBoard::SudokuBoard() {
    // Kill program if no file location provided
    fprintf(stderr, "ERROR: No file location provided.\n");
    exit(1);
}

SudokuBoard::SudokuBoard(const char *file) {
    make_board(file);
    vector_contains();
}

void SudokuBoard::make_board(const char *file) {
    // Declare buffer open file
    char buffer[10];
    FILE *in = fopen(file, "r");
    // Return error if file invalid
    if(!in) {
        fprintf(stderr, "ERROR: No file location provided.\n");
        exit(1);
    }
    // Get first entry and convert to uint
    uint16_t first = fgetc(in);
    // Clear newline
    fgetc(in);
    atoui(&first);
    // If sqrt(first) != 2/3 invalid size, report error
    if (Board::isqrt(first) == 2 || Board::isqrt(first) == 3) {
        // Instantiate new board
        auto *board = new Board(first);
        uint16_t row = 0;
        // Parse file to EOF
        while (fgets(buffer, (board->dim + board->cell), in)) {
            uint16_t V, R, C;
            // Read line
            for (int i = 0; i < board->dim; ++i) {
                // Get char and read into frame
                first = (unsigned char) buffer[i];
                V = (i/board->cell) + board->cell*(row/board->cell);
                R = row % board->cell;
                C = i % board->cell;
                if (first == '*') {
                    board->state[V].second[R][C] = first;
                    board->missing.emplace_back(row, i);
                } else {
                    atoui(&first);
                    // Populate board->state[cell][row][col]
                    board->state[V].second[R][C] = first;
                }
            }
            row++;
        }
        this->state = board;
        fclose(in);
    } else {
        fprintf(stderr, "\nERROR: Invalid board dimension.\n");
        fclose(in);
        exit(1);
    }
}

uint16_t Board::isqrt(uint16_t n) {
    uint16_t root, remainder, place;
    root = 0;
    remainder = n;
    place = 0x4000;
    while (place > remainder) place = place >> 2;
    while (place) {
        if (remainder >= root + place) {
            remainder = remainder - root - place;
            root = root + (place << 1);
        }
        root = root >> 1;
        place = place >> 2;
    }
    return root;
}

void SudokuBoard::atoui(uint16_t *out) {
    *out = *out - '0';
}

void SudokuBoard::vector_contains() const {
    unordered_set<uint16_t> *ref;
    uint16_t temp;
    // Loop over each sub-matrix
    for(uint16_t i = 0; this->state->dim > i; ++i) {
        for (uint16_t j = 0; j < this->state->cell; ++j) {
            for (uint16_t k = 0; this->state->cell > k; ++k) {
                ref = &this->state->state[i].first;
                temp = this->state->state[i].second[j][k];
                if(temp != '*' && ref->find(temp) == ref->end()) {
                    ref->insert(temp);
                } else if (temp != '*' && ref->find(temp) != ref->end()) {
                    fprintf(stderr, "ERROR: Invalid game board loaded.\n");
                    exit(1);
                }
            }
        }
    }
}