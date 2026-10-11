#include<iostream>
#include<algorithm>
#include<cstring>
#include<vector>
#include<string>
#include<sstream>
#include<map>
#include<set>
#include<cctype>

class term { //项
public:
    int a, b, c, d;

    term() {}
    term(int _a, int _b, int _c, int _d) { // 你可以用 term t(1, 2, 3, 4); 创建一个表示 x^2sin^3xcos^4x 的term
        a = _a; b = _b; c = _c; d = _d;
    }
    bool operator == (const term &obj) const { // 相等表示是同类项
        return b == obj.b && c == obj.c && d == obj.d;
    }
    bool operator != (const term &obj) const {
        return b != obj.b || c != obj.c || d != obj.d;
    }
    
    bool operator < (const term &obj) const {
        if (b != obj.b) return b > obj.b;
        if (c != obj.c) return c > obj.c;
        return d > obj.d;
    }
};

class poly { //多项式
public:
    int n;
    term *t;

    poly() {n = 0; t = NULL;}
    poly(int _n) { // 构造函数
        n = _n;
        t = new term[n];
    }
    poly(const poly &p) {
        n = p.n;
        t = new term[n];
        for (int i = 0; i < n; ++i) {
            t[i] = p.t[i];
        }
    }
    
    void simplify() { // 合并同类项 
        if (n <= 1) return;
        
        std::map<std::tuple<int,int,int>, int> term_map;
        for (int i = 0; i < n; ++i) {
            auto key = std::make_tuple(t[i].b, t[i].c, t[i].d);
            term_map[key] += t[i].a;
        }
        
        // Remove zero coefficient terms
        std::vector<term> terms;
        for (auto &pair : term_map) {
            if (pair.second != 0) {
                terms.push_back(term(pair.second, std::get<0>(pair.first), std::get<1>(pair.first), std::get<2>(pair.first)));
            }
        }
        
        // Sort terms
        std::sort(terms.begin(), terms.end());
        
        // Update array
        if (t != NULL) delete []t;
        n = terms.size();
        if (n > 0) {
            t = new term[n];
            for (int i = 0; i < n; ++i) {
                t[i] = terms[i];
            }
        } else {
            t = NULL;
        }
    }
    
    poly operator + (const poly &obj) const { // 对 poly 重载加法运算符
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) ans.t[i + n] = obj.t[i];
        ans.simplify();
        return ans;
    }
    
    poly operator - (const poly &obj) const { // 对 poly 重载减法运算符
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) {
            ans.t[i + n] = obj.t[i];
            ans.t[i + n].a *= -1;
        }
        ans.simplify();
        return ans;
    }
    
    poly operator * (const poly &obj) const { // 对 poly 重载乘法运算符
        if (n == 0 || obj.n == 0) return poly();
        
        poly ans(n * obj.n);
        int idx = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < obj.n; ++j) {
                ans.t[idx] = term(t[i].a * obj.t[j].a, 
                                 t[i].b + obj.t[j].b,
                                 t[i].c + obj.t[j].c,
                                 t[i].d + obj.t[j].d);
                idx++;
            }
        }
        ans.simplify();
        return ans;
    }
    
    poly& operator = (const poly &obj) {
        if (&obj == this) return *this;
        n = obj.n;
        if (t != NULL) delete []t;
        t = new term[n];
        for (int i = 0; i < n; ++i) {
            t[i] = obj.t[i];
        }
        return *this;
    }
    
    poly derivate() const {  // 注意,你需要以"只读模式"访问当前的poly,并返回一个求导后的poly
        if (n == 0) return poly();
        
        std::vector<term> terms;
        for (int i = 0; i < n; ++i) {
            int a = t[i].a, b = t[i].b, c = t[i].c, d = t[i].d;
            
            if (b > 0) {
                // derivative of ax^b = abx^(b-1)
                terms.push_back(term(a * b, b - 1, c, d));
            }
            
            if (c > 0) {
                // derivative of a*sin^c(x)*cos^d(x) with respect to sin part
                // using chain rule: c*a*sin^(c-1)(x)*cos(x)*cos^d(x)
                terms.push_back(term(a * c, b, c - 1, d + 1));
            }
            
            if (d > 0) {
                // derivative of a*sin^c(x)*cos^d(x) with respect to cos part
                // using chain rule: -d*a*sin^c(x)*cos^(d-1)(x)*sin(x)
                terms.push_back(term(-a * d, b, c + 1, d - 1));
            }
        }
        
        poly result(terms.size());
        for (size_t i = 0; i < terms.size(); ++i) {
            result.t[i] = terms[i];
        }
        result.simplify();
        return result;
    }
    
    
    
    ~poly() { // 析构函数,释放内存
        if (t != NULL) delete []t;
    }
};

class frac { //分式
public:
    poly p, q;

    frac() {}
    frac(int x) {
        p = poly(1);
        p.t[0] = term(x, 0, 0, 0);
        q = poly(1);
        q.t[0] = term(1, 0, 0, 0);
    }
    frac(term _p) {
        q = poly(1);
        q.t[0] = term(1, 0, 0, 0);
        p = poly(1);
        p.t[0] = _p;
    }
    frac(poly _p, poly _q) : p(_p), q(_q) {}

    frac operator + (const frac &obj) const {
        return frac(p * obj.q + q * obj.p, q * obj.q); // p1/q1+p2/q2=(p1*q2+p2*q1)/(q1*q2)
    }
    
    frac operator - (const frac &obj) const {
        return frac(p * obj.q - q * obj.p, q * obj.q); // p1/q1-p2/q2=(p1*q2-p2*q1)/(q1*q2)
    }
    
    frac operator * (const frac &obj) const {
        return frac(p * obj.p, q * obj.q); // p1/q1*p2/q2=(p1*p2)/(q1*q2)
    }
    
    frac operator / (const frac &obj) const {
        return frac(p * obj.q, q * obj.p); // p1/q1/p2/q2=(p1*q2)/(q1*p2)
    }
    
    frac derivate() const { // 对分式求导
        // (p/q)' = (p'*q - q'*p)/(q*q)
        return frac(p.derivate() * q - q.derivate() * p, q * q);
    }
    
    std::string term_to_string(const term &t) const {
        std::string result;
        
        // Handle coefficient
        if (t.b == 0 && t.c == 0 && t.d == 0) {
            // Constant term
            return std::to_string(t.a);
        }
        
        // Handle sign and coefficient
        if (std::abs(t.a) == 1) {
            if (t.a == -1) result += "-";
            // coefficient 1 is omitted for non-constant terms
        } else {
            result += std::to_string(t.a);
        }
        
        // Handle x^b
        if (t.b > 0) {
            result += "x";
            if (t.b > 1) result += "^" + std::to_string(t.b);
        }
        
        // Handle sin^c(x)
        if (t.c > 0) {
            result += "sin";
            if (t.c > 1) result += "^" + std::to_string(t.c);
            result += "x";
        }
        
        // Handle cos^d(x)
        if (t.d > 0) {
            result += "cos";
            if (t.d > 1) result += "^" + std::to_string(t.d);
            result += "x";
        }
        
        return result;
    }
    
    std::string poly_to_string(const poly &p) const {
        if (p.n == 0) return "0";
        if (p.n == 1) return term_to_string(p.t[0]);
        
        std::string result = "(";
        for (int i = 0; i < p.n; ++i) {
            if (i > 0) {
                result += (p.t[i].a > 0) ? "+" : "-";
            } else {
                if (p.t[i].a < 0) result += "-";
            }
            
            term temp = p.t[i];
            temp.a = std::abs(temp.a);
            
            result += term_to_string(temp);
        }
        result += ")";
        return result;
    }
    
    void output() {
        // Check if numerator is 0
        if (p.n == 0) {
            std::cout << "0" << std::endl;
            return;
        }
        
        // Check if denominator is 1
        if (q.n == 1 && q.t[0].a == 1 && q.t[0].b == 0 && q.t[0].c == 0 && q.t[0].d == 0) {
            std::string result = poly_to_string(p);
            // Remove outer parentheses if present and this is the main output
            if (result.length() > 2 && result[0] == '(' && result[result.length()-1] == ')') {
                result = result.substr(1, result.length() - 2);
            }
            std::cout << result << std::endl;
            return;
        }
        
        // Normal case: (numerator)/(denominator)
        std::cout << poly_to_string(p) << "/" << poly_to_string(q) << std::endl;
    }
};

// Expression parser
class ExpressionParser {
private:
    std::string s;
    int pos;
    int n;
    
    void skip_whitespace() {
        while (pos < n && isspace(s[pos])) pos++;
    }
    
    char peek() {
        skip_whitespace();
        return pos < n ? s[pos] : '\0';
    }
    
    char get() {
        skip_whitespace();
        return pos < n ? s[pos++] : '\0';
    }
    
    bool is_digit(char c) {
        return c >= '0' && c <= '9';
    }
    
    bool is_letter(char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    }
    
    int parse_number() {
        skip_whitespace();
        int result = 0;
        bool negative = false;
        bool has_digits = false;
        
        if (peek() == '-') {
            negative = true;
            get();
        }
        
        while (is_digit(peek())) {
            result = result * 10 + (get() - '0');
            has_digits = true;
        }
        
        if (!has_digits) {
            return negative ? -1 : 1;
        }
        
        return negative ? -result : result;
    }
    
    term parse_term() {
        skip_whitespace();
        int a = 1, b = 0, c = 0, d = 0;
        
        // Parse coefficient if present
        if (is_digit(peek()) || peek() == '-') {
            a = parse_number();
        }
        
        // Parse variable parts - a term can have multiple variable parts multiplied together
        while (pos < n) {
            skip_whitespace();
            if (peek() == 'x') {
                get();
                b = 1;
                if (peek() == '^') {
                    get();
                    b = parse_number();
                }
            } else if (peek() == 's') {
                // sinx or sin^nx
                get(); // s
                get(); // i
                get(); // n
                c = 1;
                if (peek() == '^') {
                    get();
                    c = parse_number();
                }
                if (peek() == 'x') {
                    get(); // x
                }
            } else if (peek() == 'c') {
                // cosx or cos^nx
                get(); // c
                get(); // o
                get(); // s
                d = 1;
                if (peek() == '^') {
                    get();
                    d = parse_number();
                }
                if (peek() == 'x') {
                    get(); // x
                }
            } else {
                break;
            }
        }
        
        return term(a, b, c, d);
    }
    
    poly parse_poly() {
        // Only parse a single term - let higher level parsers handle + and -
        term t = parse_term();
        poly result(1);
        result.t[0] = t;
        return result;
    }
    
    frac parse_factor() {
        skip_whitespace();
        
        if (peek() == '(') {
            get(); // '('
            frac result = parse_expression();
            skip_whitespace();
            if (peek() == ')') {
                get(); // ')'
            }
            return result;
        } else {
            // Parse a polynomial
            poly p = parse_poly();
            poly denominator(1);
            denominator.t[0] = term(1, 0, 0, 0);
            return frac(p, denominator);
        }
    }
    
    frac parse_term_expr() {
        frac result = parse_factor();
        
        while (pos < n) {
            skip_whitespace();
            char op = peek();
            if (op == '*') {
                get();
                frac next_factor = parse_factor();
                result = result * next_factor;
            } else if (op == '/') {
                get();
                frac next_factor = parse_factor();
                result = result / next_factor;
            } else {
                break;
            }
        }
        
        return result;
    }
    
    frac parse_expression() {
        frac result = parse_term_expr();
        
        while (pos < n) {
            skip_whitespace();
            char op = peek();
            if (op == '+') {
                get();
                frac next_term = parse_term_expr();
                result = result + next_term;
            } else if (op == '-') {
                get();
                frac next_term = parse_term_expr();
                result = result - next_term;
            } else {
                break;
            }
        }
        
        return result;
    }
    
public:
    ExpressionParser(const std::string &str) : s(str), pos(0), n(str.length()) {}
    
    frac parse() {
        return parse_expression();
    }
};

void solve(char *s, int n) {
    std::string str(s, n);
    ExpressionParser parser(str);
    frac expression = parser.parse();
    
    // Output the simplified expression
    expression.output();
    
    // Output the derivative
    frac derivative = expression.derivate();
    derivative.output();
}

int main() {
    std::string str;
    std::cin >> str;
    int n = str.length();
    char *s = new char[n + 2]{0};
    for (int i = 0; i < n; ++i) s[i] = str[i];
    solve(s, n);
    delete []s;
    return 0;
}