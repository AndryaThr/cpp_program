/*******************************************
 *                  4-LINE                 *
 *******************************************/

#include <iostream>

void stack_pile(int, int, int*);
bool is_valid(int, int*);

using namespace std;

class Line4 {
public:
    Line4();
    ~Line4();
    void init();
    void playing();
    void display_board();
    void display_game();
    int game_over(int);

private:
    int** _board;
    int _turn;
};

int main() {
    Line4 new_game;
    new_game.playing();
    return 0;
}

Line4::Line4() {
    this->_board = new int*[7];
    this->_turn = 1;
    for(int i(0); i<7; i++) {
        this->_board[i] = new int[10];
    }
    for(int i(0); i<7; i++)
        for(int j(0); j<10; j++)
            this->_board[i][j] = 0;
}

void Line4::playing() {
    int cell(0);
    while(true) {
        if(cell!=0) {
            if(is_valid(10, this->_board[cell-1])) {
                cout << ">>>> Last move: on cell " << cell << endl;
            } else {
                cout << ">>>> Invalid move, check if cell_" << cell << " is full or not" << endl;
            }
        } else {
            cout << ">>>> Player 1, it's your turn !" << endl;
        }

        this->display_game();
        cout << "Player " << this->_turn << " : ";    cin >> cell;
        stack_pile(this->_turn, 10, this->_board[cell-1]);
        if(this->game_over(this->_turn)==1) {
            cout << "Congrats, you win!!" << endl;
            break;
        } else if(this->game_over(this->_turn)==0) {
            cout << "No empty cell left!!" << endl;
        }
        if(is_valid(10, this->_board[cell-1])) {
            if(this->_turn==1)
                this->_turn = 2;
            else
                this->_turn = 1;
        }
        system("cls");
    }
}

void Line4::display_game() {
    cout << endl;

    for(int i(0); i<10; i++) {
        for(int j(0); j<7; j++) {
            cout << "+---";
        }   cout << "+";
        cout << endl;
        for(int j(0); j<7; j++) {
            cout << "| ";
            if(this->_board[j][i]==0)
                cout << " ";
            else if(this->_board[j][i]==1)
                cout << "X";
            else
                cout << "O";
            cout << " ";
        }   cout << "|";
        cout << endl;
    }

    for(int j(0); j<7; j++) {
        cout << "+---";
    }   cout << "+";
    cout << endl;
    for(int j(0); j<7; j++) {
        cout << "| " << j+1 << " ";
    }   cout << "|";    cout << endl;
    for(int j(0); j<7; j++) {
        cout << "+---";
    }   cout << "+";
    cout << endl;
}

void Line4::display_board() {
    for(int i(0); i<10; i++) {
        for(int j(0); j<7; j++) {
            cout << this->_board[j][i] << "\t";
        }
        cout << endl;
    }
}

int Line4::game_over(int player) {
    return -1;
}

Line4::~Line4() {
    for(int i(0); i<7; i++)
        delete[] this->_board[i];

    delete[] this->_board;
}

/*------------------------------------------------*/
/*------------------------------------------------*/
/*------------------------------------------------*/
/*------------------------------------------------*/
/*------------------------------------------------*/

void stack_pile(int player, int pile_size, int* pile) {
    int fmax(0);
    for(int i(0); i<pile_size; i++) {
        if(pile[i]!=0) {
            fmax = pile[i];
            pile[i-1] = player;
        }
    }
    if(fmax==0)
        pile[pile_size-1] = player;
}

bool is_valid(int pile_size, int* pile) {
/** \brief check if move is valid
 *
 * \out false : pile is full
 * \out true : valid move
 *
 */

    return (pile[0]==0);
}
