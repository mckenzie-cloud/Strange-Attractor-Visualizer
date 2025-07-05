
namespace Attractor {
    // F(t, y)
    void Lorenz(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        Ks[0] = 10.0 * (y - x);
        Ks[1] = x * (28.0 - z) - y;
        Ks[2] = x * y - (8.0/3.0) * z;
    }

    void Thomas(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        const long double b = 0.208186;
        Ks[0] = std::sinl(y) - b * x;
        Ks[1] = std::sinl(z) - b * y;
        Ks[2] = std::sinl(x) - b * z;
    }

    void Aizawa(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double a = 0.95, b = 0.7, c = 0.6, d = 3.5, e = 0.25, f = 0.1;
        long double x = Xn[0], y = Xn[1], z = Xn[2];

        Ks[0] = (z - b)*x - d*y;
        Ks[1] = d*x + (z - b)*y;
        Ks[2] = c + a*z - (std::powl(z, 3) / 3.0) - (x*x + y*y) * (1 + e*z) + f*z*std::powl(x, 3.0);
    }

    void Dadras(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks) {
        long double a = 3.0, b = 2.7, c = 1.7, d = 2.0, e = 9.0;
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        Ks[0] = y - a*x + b*y*z;
        Ks[1] = c*y - x*z + z;
        Ks[2] = d*x*y - e*z;
    }

    void Rossler(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        long double a = 0.2, b = 0.2, c = 5.7;
        Ks[0] = -y - z;
        Ks[1] = x + a * y;
        Ks[2] = b + z * (x - c);
    }

    void SprottB(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        long double a = 0.4, b = 1.2, c = 1.0;
        Ks[0] = a * y * z;
        Ks[1] = x - b * y;
        Ks[2] = c - x * y;
    }

    void Arneodo(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        long double a = -5.5, b = 3.5, c = -1.0;
        Ks[0] = y;
        Ks[1] = z;
        Ks[2] = -a * x - b * y - z + c * (x*x*x);
    }

    void Lorenz83(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks) {
        long double a = 0.95, b = 7.91, f = 4.83, g = 4.66;
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        
        Ks[0] = -a*x - y*y - z*z + a*f;
        Ks[1] = -y + x*y - b*x*z + g;
        Ks[2] = -z + b*x*y + x*z;
    }

    void Halvorsen(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        const long double a = 1.4;
        Ks[0] = -a*x- 4.0 * y - 4.0 * z - y*y;
        Ks[1] = -a*y - 4.0 * z - 4.0 * x- z*z;
        Ks[2] = -a*z - 4.0 * x- 4.0 * y - x*x;
    }

    void Sprott_LinzF(long double t, std::vector<long double> &Xn, std::vector<long double> &Ks)
    {
        long double x = Xn[0], y = Xn[1], z = Xn[2];
        long double a = 0.5;
        Ks[0] = y + z;
        Ks[1] = -x + a * y;
        Ks[2] = x * x - z;
    }
};