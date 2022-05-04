#include <iostream>
#include <ctime>
#include <cctype>
#include <vector>

using namespace std;

template<class C> void swap_element(C&, C&);

class Game {
private:
    int **board;
    int board_size;
    int moves;
    char index[4];

public:
    Game();
    ~Game();

    void display_board();
    bool is_valid(char);

    bool game_over();
    void solve();

    void init();
    void move_piece(char);
    bool validate_move(int, int, char);

    void play();
};

int main()
{
    Game newgame;
    newgame.play();
}

void Game::play() {
    init();
    char mv(0);
    while(true) {
        display_board();
        if(game_over()) {
            cout << "CONGRATS YOU DID IT!! " << endl;
            break;
        }
        while(true) {
            cout << "move ( z - q - s - d ) : "; cin >> mv;
            switch(mv) {
                case 'z':
                mv = 'u';
                break;
                case 'd':
                mv = 'r';
                break;
                case 's':
                mv = 'd';
                break;
                case 'q':
                mv = 'l';
                break;
                case 'x':
                exit(0);
            }
            if(is_valid(mv))
                break;
        }
        this->moves++;
        move_piece(mv);
        system("cls");
    }
}

Game::Game() {
    this->board_size = 3;
    this->moves = 0;
    this->board = new int* [this->board_size];
    index[0] = 'u';
    index[1] = 'l';
    index[2] = 'r';
    index[3] = 'd';

    for(int i(0); i<this->board_size; i++) {
        this->board[i] = new int[this->board_size];
        for(int j(0); j<this->board_size; j++) {
            this->board[i][j] = i*3 + j + 1;
        }
    }
}

Game::~Game() {
    for(int i(0); i<this->board_size; i++)
        delete[] this->board[i];

    delete[] this->board;
}

void Game::init() {
    int shuffle = 300;
    int mv = 0;
    srand(time(0));
    for(int i(0); i<shuffle; ++i) {
        while(true) {
            mv = rand() % 4;
            if(is_valid(index[mv]))
                break;
        }
        this->move_piece(index[mv]);
    }
}

void Game::display_board() {
    cout << "moves : " << this->moves << endl;
    for(int i(0); i<this->board_size; i++) {
        cout << "++-----+-----+-----++" << endl;
        cout << "||";
        for(int j(0); j<this->board_size; j++) {
            cout << "  ";
            (board[i][j] != 9) ? cout << board[i][j] : cout << " ";
            cout << "  |";
        }
        cout << "|" << endl;
    }
    cout << "++-----+-----+-----++" << endl;
}

void Game::move_piece(char direction) {
    int i(0), j(0);
    bool end_loop(false);
    for(i=0; i<board_size; i++) {
        for(j=0; j<board_size; j++) {
            if(board[i][j] == 9) {
                end_loop = true;
                break;
            }
        }
        if(end_loop)
            break;
    }

    if(!validate_move(i,j,direction)) {
        cout << "INVALID MOVE " << endl;
    } else {
        switch(tolower(direction)) {
            case 'd':
                swap_element(board[i][j], board[i-1][j]);
                break;
            case 'r':
                swap_element(board[i][j], board[i][j-1]);
                break;
            case 'l':
                swap_element(board[i][j], board[i][j+1]);
                break;
            case 'u':
                swap_element(board[i][j], board[i+1][j]);
                break;
        }
    }
}

void Game::solve() {
    vector<vector<int>> neighbors;
    for(int i(0); i<9; i++) {
        int l=0;
        if(i%2==0 && i!=4)
            l=2;
        else if(i%2==1)
            l=3;
        else
            l=4;
        neighbors.push_back(vector<int>(l));
    }
}

bool Game::is_valid(char dep) {
    int i(0), j(0);
    bool end_loop(false);
    for(i=0; i<board_size; i++) {
        for(j=0; j<board_size; j++) {
            if(board[i][j] == 9) {
                end_loop = true;
                break;
            }
        }
        if(end_loop)
            break;
    }
    return validate_move(i,j,dep);
}

bool Game::validate_move(int i, int j, char dep) {
    switch(tolower(dep)) {
        case 'd':
            if(i==0)
                return false;
            break;
        case 'r':
            if(j==0)
                return false;
            break;
        case 'l':
            if(j==2)
                return false;
            break;
        case 'u':
            if(i==2)
                return false;
            break;
    }
    return true;
}

bool Game::game_over() {
    for(int i(0); i<this->board_size; i++) {
        for(int j(0); j<this->board_size; j++) {
            if(board[i][j] != i*3 + j + 1) {
                return false;
            }
        }
    }
    return true;
}

template<class C>
void swap_element(C& a, C& b) {
    int c = a;
    a = b;
    b = c;
}
