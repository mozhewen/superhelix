#include <cmath>
#include <fstream>

using namespace std;


const double dx = 1E-5;
const double pi = M_PI;

class v3 {
public:
    double x, y, z;
    v3(): x(0.),y(0.),z(0.) {}
    v3(double _x, double _y, double _z): x(_x),y(_y),z(_z) {}
    inline v3 operator+(const v3 &b) const {return v3(x+b.x, y+b.y, z+b.z);}
    inline v3 operator-(const v3 &b) const {return v3(x-b.x, y-b.y, z-b.z);}
    inline v3 operator*(const double f) const {return v3(x*f, y*f, z*f);}
    inline v3 operator/(const double f) const {return v3(x/f, y/f, z/f);}
    inline double norm() const {return sqrt(x*x+y*y+z*z);}
    inline v3 & normalize() {
        double norm = sqrt(x*x+y*y+z*z);
        x /= norm; y /= norm; z /= norm;
        return *this;
    }
    friend inline v3 cross(const v3 &a, const v3 &b) {
        return v3(
            a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x
        );
    }
};


v3 d(v3 f(double), double x)
{
    return (f(x+dx) - f(x-dx)) / (2.*dx);
}

v3 d2(v3 f(double), double x)
{
    return (f(x+dx) - f(x)*2. + f(x-dx)) / (dx*dx);
}

inline double sqr(double x) {
    return x*x;
}

v3 helixOn(v3 f(double), double x, double A, double omega, v3 n = v3())
{
    v3 t = d(f, x).normalize();
    if (n.norm() < 1E-6)
        n = d2(f, x).normalize();
    v3 n1 = cross(t, n).normalize();
    return f(x) + (n*cos(omega*x) + n1*sin(omega*x))*A;
}


const int N = 100000;
v3 p[N];
int i;
double xi;

int main()
{
    ofstream fout("out.txt");

    for (int i = 0; i < N; i++) {
        xi = (double)i/(double)N;
        p[i] = helixOn(
            [](double x) {return helixOn(
                [](double x) {return helixOn(
                    [](double x) {return v3(0,0,1)*4*x;}, x,
                    1, 2*pi*4,                      // A, omega
                    v3(1,0,0)
                );}, x,
                4*sqrt(1+sqr(2*pi))/160, 2*pi*160   // A, omega
            );}, xi,
            4*(1+sqr(2*pi))/6400, 2*pi*6400         // A, omega: Only approximately correct.
        );
    }

    for (int i = 0.; i < N; i++)
        fout << p[i].x << ' ' << p[i].y << ' ' << p[i].z << endl;

    fout.close();

    return 0;
}
