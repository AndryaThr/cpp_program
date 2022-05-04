#include <iostream>
#include <string>
#include <cctype>

using namespace std;

void compute(string);

int main()
{
    string expression("");

    cout << "Entrez l'expression a calculer : ";    getline(cin, expression);
    cout << expression << endl;
    compute(expression);

    return 0;
}

void compute(string expression) {
    string expr1, expr2;
    char sign('#');
    for(int i(0); i<expression.length(); i++) {

        if(isdigit(expression[i])) {
            if(sign=='#') {
                expr1.push_back(expression[i]);
            } else {
                expr2.push_back(expression[i]);
            }
            expression[i] = ' ';
        } else if(!isdigit(expression[i]) && !isspace(expression[i])) {
            sign = expression[i];
        }

    }
    float res(0);
    switch(sign) {
    case '+':
        res = stof(expr1) + stof(expr2);
        break;
    case '*':
        res = stof(expr1) * stof(expr2);
        break;
    case '-':
        res = stof(expr1) - stof(expr2);
        break;
    case '/':
        res = stof(expr1) / stof(expr2);
        break;
    default:
        cout << "Unknown sign" << endl;
        exit(-1);
    }

    cout << expr1 << " " << sign << " " << expr2 << " = " << res << endl;
}

