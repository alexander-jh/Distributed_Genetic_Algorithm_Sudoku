// SudokuBoard.hpp
//
// Header files for the game board representation.
//
// Version:    C++11
// Date:       11/19/2020
// Author:     Alex Hoke
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <unordered_set>
#include <queue>
#include <cstdlib>
#include <random>
#include <unordered_map>

using namespace std;

// Struct for the sudoku board
//  @param  state   -   vector containing dim elements of cell x cell
//                      sub-matrices
//  @param  missing -   vector of reference to missing positions
//  @param  dim     -   overall dimension of the sudoku board
//  @param  cell    -   dimension of of sub matrices, sqrt(dim)
typedef struct Board {
    vector<pair<unordered_set<uint16_t>, vector<vector<uint16_t>>>> state;
    vector<pair<uint16_t,uint16_t> *>                               missing;
    uint16_t                                                        dim;
    uint16_t                                                        cell;
    // Constructor to ensure proper vector space is allocated
    // @param   dim -   Overall dimension of the sudoku board
    explicit Board(uint16_t dim);
    ~Board();
    // Unsigned integer square root function
    // @param   n   -   16-bit unsigned integer
    static uint16_t isqrt(uint16_t n);
} Board;

// Representation of the initial game board to be loaded into the
// optimization algorithms. Verifies that initial input is solvable
class SudokuBoard {
public:
    // Constructor for SudokuBoard when appropriate input is supplied
    // @param file  -   location to a file
    SudokuBoard(const char *file);
    // Constructor which prompts user for file input if not supplied
    SudokuBoard();
    // Destructor
    inline ~SudokuBoard(){delete state;};
    // Current board state broken down into a vector of 4 2x2 sub-matrices
    // to make verification of the initial board space easier
    Board *state{};
private:
    // Parses file input following verification of the first character
    // @param *file -   open file pointer to parse
    void make_board(const char *file);
    // Converts character to 16 bit unsigned integer
    // @param *out  -   character to convert
    static void atoui(uint16_t *out);
    // Fill set with elements observed in each vector to verify initial state
    void vector_contains();
};