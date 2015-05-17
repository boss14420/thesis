#ifndef __COMMON_HPP__
#define __COMMON_HPP__

#include <fstream>
#include <cmath>

template <typename T>
void printColor(std::ostream &os, T value)
{
    unsigned char color = std::floor(255-255*value);
//    std::uint8_t rgb[] = { color, color, color };
    os.write(reinterpret_cast<char const*>(&color), 1);
//    os << std::floor(255-255*value) << ' '
//       << std::floor(255-255*value) << ' '
//       << std::floor(255-255*value) << ' ';
}

template <typename T>
void exportPixmap(T const *c, std::size_t width, std::size_t height, std::size_t stride, char const *filename)
{
    std::ofstream ofs(filename); //, std::ios::binary);
    ofs << "P5\n";
    ofs << width + 1 << ' ' << height + 1 << "\n255\n";

    for (int j = height; j >= 0; --j) {
        for (std::size_t i = 0; i <= width; ++i) {
            printColor(ofs, std::fabs(c[i*stride + j]));
        }
    }

//    T const *prow = c + width * (height - 1);
//    for (; prow != c-width; prow -= width) {
//        for (std::size_t x = 0; x != width; ++x)
//            printColor(ofs, prow[x]);
//    }
    ofs.close();
}

template <typename T>
void exportMatlab(T const *c, std::size_t width, std::size_t height, std::size_t stride, char const *filename)
{
    std::ofstream ofs(filename);
//    ofs << "# name: u\n"
//        << "# type: matrix\n"
//        << "# rows: " << HEIGHT+1 << '\n'
//        << "# columns: " << WIDTH+1 << '\n';
    for (std::size_t j = 0; j <= height; ++j) {
        for (std::size_t i = 0; i <= width; ++i) {
            ofs << c[i*stride+j] << " ";
        }
        ofs << '\n';
    }
}

//template <typename T>
//void exportValue(T const *v)
//{
//    for (int j = HEIGHT; j >= 0; --j) {
//        for (int i = 0; i <= WIDTH; ++i) {
//            std::cout << std::setw(5) << v[j*STRIDE+i] << ' ';
//        }
//        std::cout << '\n';
//    }
//}

template <typename T>
T stddev(T const *c1, T const *c2, std::size_t width, std::size_t height) {
    T sum = 0;
    for (T const *c1tmp = c1; c1tmp != c1 + width*height; ++c1tmp, ++c2) {
        sum += (*c1tmp - *c2) * (*c1tmp - *c2);
    }
    return std::sqrt(sum / (width * height));
}

#endif
