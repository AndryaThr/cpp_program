/****************************************************
 *                  TOUR D'HANOI                    *
 ****************************************************/

#include <iostream>
#include <cstdlib>

#include<chrono>
#include<thread>

void move_to(int, int*, int*);
int is_valid(int, int*, int*);
void display_table(int, int*);

using namespace std;

class Hanoi_tower {
public:
    Hanoi_tower();
    Hanoi_tower(int);
    ~Hanoi_tower();
    void init();
    void display_game();
    void display_board();
    void playing();
    bool game_over();
    void solve(int);
    void solve_2(int, int, int*, int*, int*);

private:
    int _N;
    int** _board;
};

int main() {
    cout << "TOWER OF HANOI\n" << endl;

    int level(3);
    while(true) {
        cout << "Choose your level (3 ~ 7) : ";
        cin >> level;
        if(level >= 2 && level <= 7)
            break;
    }

    // create game
    Hanoi_tower HGame(level);
    HGame.init();
    //HGame.playing();
    HGame.solve(level);

    return 0;
}

Hanoi_tower::Hanoi_tower() {
    this->_N = 3;
}

Hanoi_tower::Hanoi_tower(int level) {
    this->_N = level;
}

void Hanoi_tower::init() {
    this->_board = new int*[3];
    for(int i(0); i<3; i++) {
        this->_board[i] = new int[this->_N];
    }

    for(int i(0); i<this->_N; i++) {
        this->_board[0][i] = i+1;
        this->_board[1][i] = 0;
        this->_board[2][i] = 0;
    }
}

void Hanoi_tower::display_game() {
    cout << endl;
    for(int j(0); j<3; j++)
            cout << "\t|";

    cout << endl;
    for(int i(0); i<this->_N; i++) {
        for(int j(0); j<3; j++) {
            cout << "\t";
            if(this->_board[j][i]==0)
                cout << "|";
            else
                cout << this->_board[j][i];
        }
        cout << endl;
    }
    for(int i(0); i<33; i++)
        cout << "-";
    cout << endl;
}

void Hanoi_tower::display_board() {
    // cout << "This is the board table : " << endl;
    for(int i(0); i<3; i++) {
        for(int j(0); j<this->_N; j++)
            cout << this->_board[i][j] << "\t";
        cout << endl;
    }
}

void Hanoi_tower::playing() {
    int from(0),
        to(0),
        attempt(0),
        validation(0);

    while(true) {
        system("cls");
        cout << "TOWER OF HANOI\n" << endl;

        if(from==0 && to==0) {
            cout << ">>>> Starting game ..." << endl;
        } else {
            if(validation==-1)
                cout << ">>>> Nothing to move on tower n_" << from << endl;
            else if(validation==0)
                cout << ">>>> Invalid move, check if the tower n_" << to << " is full or have a larger number on top !" << endl;
            else if(validation==1)
                cout << ">>>> last move: tower n_" << from << " to tower n_" << to << endl;
        }

        cout << "attempt : " << attempt << endl;
        this->display_game();

        while(true) {
            cout << "From (1 ~ 3): ";  cin >> from;
            cout << "To (1 ~ 3): ";  cin >> to;
            if((from <= 0 || from >3) || (to <= 0 || to >3) || from == to) {
                cout << ">>>> Retry !!" << endl;
            } else {
                cout << endl;
                break;
            }
        }

        validation = is_valid(this->_N, this->_board[from-1], this->_board[to-1]);
        move_to(this->_N, this->_board[from-1], this->_board[to-1]);
        attempt++;

        if(this->game_over()) {
            cout << "Congrats!! You win on " << attempt << " moves !!" << endl;
            break;
        }
    }
}

bool Hanoi_tower::game_over() {
    for(int i(0); i<this->_N; i++) {
        if(this->_board[0][i]!=0 || this->_board[1][i]!=0)
            return false;
    }
    return true;
}

Hanoi_tower::~Hanoi_tower () {
    for(int i(0); i<3; i++) {
        delete this->_board[i];
    }
    delete this->_board;
}

void Hanoi_tower::solve_2(int N, int level, int* start, int* mid, int* stop) {
    if(game_over()) {
        cout << "Solved !!" << endl;
        return;
    } else {
        if(level == 0) {
            return;
        }
        solve_2(N, level-1, start, stop, mid);
        display_game();
        this_thread::sleep_for(chrono::milliseconds(700));
        system("cls");
        move_to(N, start, stop);
        solve_2(N, level-1, mid, start, stop);
    }
}

void Hanoi_tower::solve(int level) {
    display_game();
    solve_2(level, _N, this->_board[0], this->_board[1], this->_board[2]);
}

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
/* -------------------------------------------------------- */

void move_to(int N, int* from, int* to) {
    int ftop(0),
        fid(0),
        ttop(0);
    for(int i(0); i<N; i++) {
        if(from[i]!=0) {
            ftop = from[i];
            fid = i;
            break;
        }
    }
    if(ftop==0) {
        cout << "Invalid move !" << endl;
        return;
    }

    for(int j(0); j<N; j++) {
        if(to[j]!=0) {
            ttop = to[j];
            if(j==0 || ttop < ftop) {
                cout << "Invalid move !" << endl;
                return;
            } else {
                to[j-1] = ftop;
                from[fid] = 0;
            }
            break;
        } else if(j==N-1 && to[j]==0) {
            to[j] = ftop;
            from[fid] = 0;
        }
    }
}

int is_valid(int N, int* from, int* to) {

/** \brief check if move is valid
 *
 * \out -1: from is empty
 * \out 0 : ftop is larger than ttop
 * \out 1 : valid move
 *
 */

    int ftop(0),
        ttop(0);

    for(int i(0); i<N; i++) {
        if(from[i]!=0) {
            ftop = from[i];
            break;
        }
    }
    if(ftop==0) {
        return -1;
    }

    for(int j(0); j<N; j++) {
        if(to[j]!=0) {
            ttop = to[j];
            if(j==0 || ttop < ftop) {
                return 0;
            }
            break;
        }
    }
    return 1;
}

void display_table(int n, int* tab) {
    for(int i(0); i<n; i++) {
        cout << tab[i] << "\t";
    }
    cout << endl;
}
