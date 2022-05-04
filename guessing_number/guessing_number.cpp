#include <iostream>

using namespace std;

void guess(int, int);

int main()
{
    guess(1, 1024);
    return 0;
}

void guess(int mn, int mx) {
    if(mx-mn<=1) {
        cout << "le nombre est : " << mn << endl;
        return;
    }

    char choice(' ');
    int half = (mn + mx)/2;
    cout << "greater or equal to " << half << " (y or n) : "; cin >> choice;
    if(choice == 'y') {
        guess(half, mx);
    } else {
        guess(mn, half);
    }
}
