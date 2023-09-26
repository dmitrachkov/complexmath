#ifndef __complex_math
#define __complex_math

/*
 * ATTENTION!
 * 
 * These are functions to operate with complex numbers.
 * I wrote this to use in Houdini.
 * Easy to read. Easy to modify for your needs.
 *
 * YOU ARE ALLOWED TO:
 * Copy, extend, convert, share, fix, use.
 *
 * THANK YOU FOR YOUR ATTENTION!
*/


//https://en.wikipedia.org/wiki/Complex_number


// Modulus of a complex number (this is absolute value
// or magnitude or length of the vector)
float cmod(vector2 z)
{
    // z = a + bi
    // |z| = sqrt(a^2 + b^2)
    // sqrt(z[0] * z[0] + z[1] * z[1]);
    return length(z);
}

// Argumnet of a complex number (phase or angle from x-axis)
float carg(vector2 z)
{
    return atan2(z[1], z[0]);
}

// Convert complex number to polar coordinates
vector2 ctop(vector2 z)
{
    float m = length(z);
    if (m == 0) return set(0, 0);
    return set(m, atan2(z[1], z[0]));
}

// Convert polar coordinates to complex number
vector2 ptoc(vector2 p)
{   
    // z = a + bi
    // a = r * cos(θ); b = r * sin(θ);
    float re = p[0] * cos(p[1]);
    float im = p[0] * sin(p[1]);
    return set(re, im);
}

// Addition of complex numbers
vector2 cadd(vector2 z; vector2 w)
{
    return set(z[0] + w[0], z[1] + w[1]);
}

// Subtraction of complex numbers
vector2 csub(vector2 z; vector2 w)
{
    return set(z[0] - w[0], z[1] - w[1]);
}

// Multiplication of complex numbers
vector2 cmult(vector2 z; vector2 w)
{
    // z = a + bi; w = c + di;
    // z * w = (ac - bd) + (ad + bc)i
    float re = z[0] * w[0] - z[1] * w[1];
    float im = z[0] * w[1] + z[1] * w[0];
    return set(re, im);
}

// Multiplication of a complex number by an imaginary number
vector2 cmulti(vector2 z; float c)
{
    // z = a + bi
    // z * i = (a + bi) = a*i + b*i^2 = ai - b = -b + ai
    float re = -z[1] * c;
    float im = z[0] * c;
    return set(re, im);
}

// Division of complex numbers
vector2 cdiv(vector2 z; vector2 w)
{
    // z = a + bi; w = c + di;
    // z / w = (z * w) / (w * w)
    // = ((ac + bd) / (c^2 + d^2)) + ((bc - ad) / (c^2 + d^2))i
    float d = w[0] * w[0] + w[1] * w[1];
    if (d == 0) return set(0, 0);
    float re = (z[0] * w[0] + z[1] * w[1]) / d;
    float im = (z[1] * w[0] - z[0] * w[1]) / d;
    return set(re, im);
}

// Natural logarithm of a complex number
// https://en.wikipedia.org/wiki/Complex_logarithm
vector2 clog(vector2 z; int k)
{
    // z = r * (cos(θ) + i*sin(θ)); r = |z|
    // ln(z) = ln(r) + i(θ + 2πk)
    float r = length(z);
    if (r == 0) return set(0, 0);
    return set(log(r), atan2(z[1], z[0]) + (2 * M_PI * k));
}

// Exponentiation of a complex number to a power of a natural number
// https://en.wikipedia.org/wiki/De_Moivre%27s_formula
vector2 cnpow(vector2 z; float n)
{
    // z^n = r^n * (cos(nθ) + i*sin(nθ)); r = |z|
    float r = pow(length(z), n);
    float theta = n * atan2(z[1], z[0]);
    return set(r * cos(theta), r * sin(theta));
}

// Exponentiation of a natural number to a power of a complex number
vector2 ncpow(float n; vector2 z)
{
    // n^z = e^(z * ln(n))
    // or
    // z = a + bi
    // n^(a + bi) = n^a * n^(bi) 
    // https://en.wikipedia.org/wiki/Euler%27s_formula
    // if n = e^ln(n) and e^(bi) = cos(b) + i*sin(b) => n^(bi) = e^(i*b*ln(n)) = cos(b*ln(n)) + i*sin(b*ln(n))
    // n^a * n^(bi) = n^a * (cos(b*ln(n)) + i*sin(b*ln(n)))
    
    // n^a
    float r = pow(n, z[0]);
    // b * ln(n)
    float theta = z[1] * log(n);
    return set(r * cos(theta), r * sin(theta));
}

// Exponentiation of a complex number to a power of a complex number
vector2 ccpow(vector2 z; vector2 w; int k)
{
    // z^w = e^(w * ln(z))
    vector2 exponent = cmult(w, clog(z, k));
    return ncpow(M_E, exponent);
}

// Complex Square Root
vector2 csqrt(vector2 z)
{
    float r = length(z);
    if (r == 0) return set(0, 0);
    float theta = 0.5 * atan2(z[1], z[0]);
    return sqrt(r) * set(cos(theta), sin(theta));
}

// nth root of a complex number
// variable k=={0..n-1}
// roots for each k gives a perfect n-gon
// https://en.wikipedia.org/wiki/De_Moivre%27s_formula
vector2 cnroot(vector2 z; int n; int k)
{
    // z = r * (cos(θ) + i*sin(θ)); r = |z|
    // root = r^(1/n) * (cos((θ + 2πk) / n) + i*sin((θ + 2πk) / n))
    float r = length(z);
    if (r == 0) return set(0, 0);
    float theta = atan2(z[1], z[0]);

    float r_root = pow(r, 1.0 / n);
    float theta_root = (theta / n) + ((2.0 * k * M_PI) / n);

    float re = r_root * cos(theta_root);
    float im = r_root * sin(theta_root);

    return set(re, im);
}

#endif