#ifndef MATRIX__MATRIX_H_
#define MATRIX__MATRIX_H_
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <sstream>
#include <string>


class BigInteger {
  int sign = 1;
  std::vector<long long> arr;
  static const int BASE = 1000000;
  static const int BASE_CNT = 6;

  static int getIntLen(long long a) {
    if (a == 0) {
      return 1;
    }
    int cnt = 0;
    while (a != 0) {
      ++cnt;
      a /= 10;
    }

    return cnt;
  }

  [[nodiscard]] int getSize() const {
    return static_cast<int>(arr.size());
  }

  void normalise() {
    int new_size = getSize();
    for (int i = getSize() - 1; i >= 0; --i) {
      if (arr[i] == 0) {
        new_size = i;
      } else {
        break;
      }
    }
    arr.resize(new_size);
  }

  void rShift() {
    if (arr.empty()) {
      arr.push_back(0);
      return;
    }
    std::vector<long long> new_arr(getSize() + 1);
    std::copy(arr.begin(), arr.end(), new_arr.begin() + 1);
    new_arr[0] = 0;
    arr = std::move(new_arr);
  }

  [[nodiscard]] int cmp_wo_sign(const BigInteger& other) const {
    if (getSize() == 0 && other.getSize() == 0) {
      return 0;
    }

    if (getSize() < other.getSize()) {
      return -1;
    } else if (getSize() > other.getSize()) {
      return 1;
    } else {
      for (int i = getSize() - 1; i >= 0; --i) {
        if (arr[i] < other.arr[i]) {
          return -1;
        } else if (arr[i] > other.arr[i]) {
          return 1;
        }
      }
    }

    return 0;
  }

 public:
  [[nodiscard]] std::string toString() const {
    std::ostringstream oss;

    if (cmp_wo_sign(0) == 0) {
      oss << 0;
      return oss.str();
    }

    if (sign == -1) {
      oss << '-';
    }
    bool zero_flag = true;
    for (int i = getSize() - 1; i >= 0; --i) {
      if (zero_flag && arr[i] == 0) {
        continue;
      }
      if (!zero_flag) {
        int int_len = getIntLen(arr[i]);
        for (int j = 0; j < (BASE_CNT - int_len); ++j) {
          oss << 0;
        }
      }
      oss << arr[i];
      zero_flag = false;
    }

    return oss.str();
  }

  BigInteger() = default;

  BigInteger(int other) {
    if (other == 0) {
      return;
    }
    arr.push_back(abs(other));
    sign = (other >= 0) ? 1 : -1;
  }

  explicit BigInteger(const std::string& str) {
    if (str == "0" || str == "-0") {
      return;
    }

    long long cell = 0;
    int cnt = 0;

    for (int i = static_cast<int>(str.size()) - 1; i > 0; --i) {
      cell = (str[i] - '0') * static_cast<long long>(pow(10, cnt)) + cell;
      ++cnt;

      if (cnt == BASE_CNT) {
        arr.push_back(cell);
        cell = 0;
        cnt = 0;
      }
    }
    if (str[0] == '-') {
      sign = -1;
    } else {
      cell = (str[0] - '0') * static_cast<long long>(pow(10, cnt)) + cell;
      ++cnt;
    }

    if (cnt != 0) {
      arr.push_back(cell);
    }
  }

  std::vector<long long>::iterator data() {
    return arr.begin();
  }

  [[nodiscard]] std::vector<long long>::const_iterator data() const {
    return arr.begin();
  }

  int& getSign() {
    return sign;
  }

  const int& getSign() const {
    return sign;
  }

  explicit operator bool() const {
    return !(*this == 0);
  }

  friend bool operator==(const BigInteger& left, const BigInteger& right) {
    if (left.sign != right.sign) {
      return false;
    }
    if (left.cmp_wo_sign(right) == 0) {
      return true;
    }

    return false;
  }

  friend bool operator<(const BigInteger& left, const BigInteger& right) {
    if (left.sign != right.sign) {
      return (left.sign == -1);
    }
    if (left.cmp_wo_sign(right) == -1) {
      return (left.sign == 1);
    }

    return false;
  }

  BigInteger& operator+=(const BigInteger& other) {
    size_t max_size = std::max(arr.size(), other.arr.size());
    if (sign != other.sign) {
      sign = other.sign;
      (*this) -= other;
      sign = (sign == 1) ? -1 : 1;
      if (cmp_wo_sign(0) == 0) {
        sign = 1;
      }
    } else {
      arr.resize(max_size + 10);
      for (int i = 0; i < getSize(); ++i) {
        arr[i] += (i < other.getSize()) ? other.arr[i] : 0;
      }

      for (int i = 0; i < getSize() - 1; ++i) {
        if (arr[i] >= BASE) {
          arr[i] -= BASE;
          arr[i + 1]++;
        }
      }
      normalise();
    }
    if (cmp_wo_sign(0) == 0) {
      sign = 1;
    }
    return *this;
  }

  BigInteger& operator-=(const BigInteger& other) {
    std::vector<long long> other_arr = other.arr;

    if (sign != other.sign) {
      sign = other.sign;
      (*this) += other;
      sign = (sign == 1) ? -1 : 1;
      if (cmp_wo_sign(0) == 0) {
        sign = 1;
      }
    } else {
      if (cmp_wo_sign(other) == 0) {
        arr.resize(0);
        sign = 1;
        return *this;
      } else if (cmp_wo_sign(other) == -1) {
        other_arr = arr;
        arr = other.arr;
        sign = (sign == 1) ? -1 : 1;
      }
      size_t max_size = std::max(arr.size(), other_arr.size());
      arr.resize(max_size + 10);
      for (int i = 0; i < getSize(); ++i) {
        arr[i] -= (i < static_cast<int>(other_arr.size())) ? other_arr[i] : 0;
      }

      for (int i = 0; i < getSize() - 1; ++i) {
        if (arr[i] < 0) {
          arr[i] += BASE;
          arr[i + 1]--;
        }
      }

      normalise();
      if (cmp_wo_sign(0) == 0) {
        sign = 1;
      }
    }

    return *this;
  }

  BigInteger& operator*=(const BigInteger& other);

  BigInteger& operator/=(const BigInteger& other);

  BigInteger& operator%=(const BigInteger& other);

  BigInteger& operator/=(int other) {
    long long carry = 0;
    for (int i = getSize() - 1; i >= 0; --i) {
      long long cur = arr[i] + carry * 1ll * BASE;
      arr[i] = cur / other;
      carry = cur % other;
    }
    normalise();

    return *this;
  }

  BigInteger operator*(const BigInteger& right) const {
    BigInteger res;

    if (cmp_wo_sign(0) == 0 || right.cmp_wo_sign(0) == 0) {
      return res;
    }

    if (sign != right.sign) {
      res.sign = -1;
    }

    res.arr.resize(getSize() + right.getSize() + 1);
    for (int i = 0; i < getSize(); i++) {
      for (int j = 0; j < right.getSize(); j++) {
        res.arr[i + j] += (arr[i] * right.arr[j]);
      }
    }

    for (int i = 0; i < res.getSize() - 1; ++i) {
      long long p_to_next = res.arr[i] / BASE;
      if (p_to_next != 0) {
        res.arr[i] %= BASE;
        res.arr[i + 1] += p_to_next;
      }
    }
    res.normalise();

    return res;
  }

  friend BigInteger operator/(const BigInteger& left, const BigInteger& right) {
    BigInteger right_c = right;
    right_c.sign = 1;
    BigInteger result, cursor;
    result.arr.resize(left.arr.size());
    for (int i = left.getSize() - 1; i >= 0; --i) {
      cursor.rShift();
      cursor.arr[0] = left.arr[i];
      cursor.normalise();

      int x_value = 0, b_left = 0, b_right = BASE;
      while (b_left <= b_right) {
        int middle = (b_left + b_right) / 2;
        BigInteger temp = right_c * middle;

        if (temp < cursor || temp == cursor) {
          x_value = middle;
          b_left = middle + 1;
        } else {
          b_right = middle - 1;
        }
      }

      result.arr[i] = x_value;
      cursor -= right_c * x_value;
    }

    result.sign = (left.sign != right.sign) ? -1 : 1;
    result.normalise();
    return result;
  }

  BigInteger& operator++() {
    *this += 1;
    return *this;
  }

  BigInteger operator++(int);

  BigInteger& operator--() {
    *this -= 1;
    return *this;
  }

  BigInteger operator--(int);

  BigInteger operator-() const {
    BigInteger res = *this;
    res.sign = (res.sign == 1) ? -1 : 1;

    return res;
  }
};
std::istream& operator>>(std::istream& in, BigInteger& num) {
  std::string str;
  in >> str;

  num = BigInteger(str);

  return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& num) {
  std::string str = num.toString();

  out << str;
  return out;
}

BigInteger operator+(const BigInteger& left, const BigInteger& right) {
  BigInteger res = left;

  res += right;

  return res;
}

BigInteger operator-(const BigInteger& left, const BigInteger& right) {
  BigInteger res = left;

  res -= right;

  return res;
}

BigInteger operator%(const BigInteger& left, const BigInteger& right) {
  return left - right * (left / right);
}

BigInteger& BigInteger::operator*=(const BigInteger& other) {
  *this = (*this) * other;

  return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& other) {
  *this = (*this) / other;

  return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& other) {
  *this = (*this) % other;

  return *this;
}

BigInteger BigInteger::operator++(int) {
  *this += 1;
  return (*this - 1);
}

BigInteger BigInteger::operator--(int) {
  *this -= 1;
  return (*this + 1);
}

bool operator!=(const BigInteger& left, const BigInteger& right) {
  return !(left == right);
}

bool operator<=(const BigInteger& left, const BigInteger& right) {
  return (left < right || left == right);
}

bool operator>(const BigInteger& left, const BigInteger& right) {
  return !(left <= right);
}

bool operator>=(const BigInteger& left, const BigInteger& right) {
  return !(left < right);
}

BigInteger operator "" _bi(const char* ch_arr, unsigned long) {
  BigInteger obj(ch_arr);

  return obj;
}

class Rational {
  static const int kCastToDoublePrecision = 53;

  BigInteger num = 0;
  BigInteger denom = 1;

  void reduct() {
    if (num == 0) {
      denom = 1;
      return;
    }
    if (num == 1) {
      return;
    }

    if (num == denom) {
      num = 1;
      denom = 1;
      return;
    }

    if (denom > num && denom % num == 0) {
      denom /= num;
      num = 1;
      return;
    }

    if (num > denom && num % denom == 0) {
      num /= denom;
      denom = 1;
      return;
    }

    BigInteger dlt = gcd(num, denom);
    if (dlt > 1) {
      num /= dlt;
      denom /= dlt;
    }
  }

  static BigInteger gcd(BigInteger a, BigInteger b) {
    if (a == 0 || b == 0) {
      return 0;
    }
    a.getSign() = 1;
    b.getSign() = 1;

    BigInteger ans = 1;
    while (a != b) {
      bool ab = a.data()[0] % 2 == 0, bb = b.data()[0] % 2 == 0;
      if (ab && bb) {
        a /= 2;
        b /= 2;
        ans *= 2;
      } else if (ab) {
        a /= 2;
      } else if (bb) {
        b /= 2;
      } else {
        if (a > b) {
          a -= b;
        } else {
          b -= a;
        }
      }
    }
    return ans * a;
  }

 public:
  [[nodiscard]] std::string toString() const {
    std::ostringstream oss;

    oss << num;
    if (denom != 1) {
      oss << '/' << denom;
    }
    return oss.str();
  }

  [[nodiscard]] std::string asDecimal(size_t precision = 0) const {
    if (precision == 0) {
      return (num / denom).toString();
    }
    BigInteger copy_of_num = num;
    std::vector<int> digits;
    BigInteger remainder = (copy_of_num / denom);
    copy_of_num.getSign() = 1;
    copy_of_num %= denom;
    int i = 0;
    ++precision;
    while (i < static_cast<int>(precision)) {
      if (copy_of_num == 0) {
        digits.push_back(0);
      } else {
        copy_of_num *= 10;
        while (copy_of_num < denom) {
          copy_of_num *= 10;
          digits.push_back(0);
        }
        if (i < static_cast<int>(precision)) {
          digits.push_back(static_cast<int>(*(copy_of_num / denom).data()));
        }
        copy_of_num %= denom;
      }
      ++i;
    }

    if (digits[precision - 1] >= 5) {
      ++digits[precision - 2];
      for (size_t j = digits.size() - 2; j > 0; --j) {
        if (digits[j] >= 10) {
          digits[j] -= 10;
          ++digits[j - 1];
        } else {
          break;
        }
      }
      if (digits[0] >= 10) {
        digits[0] -= 10;
        ++remainder;
      }

    }
    std::string out_string = (num.getSign() == -1) ? "-" : "";
    out_string += remainder.toString();
    out_string += '.';
    for (size_t j = 0; j < precision - 1; ++j) {
      out_string += std::to_string(digits[j]);
    }

    return out_string;
  }

  Rational() = default;

  Rational(const BigInteger& num_p, const BigInteger& denom_p) : num(num_p), denom(denom_p) {
    reduct();
  }

  Rational(const BigInteger& other) {
    num = other;
  }

  Rational(const int other) {
    num = other;
  }

  friend bool operator==(const Rational& left, const Rational& right) {
    return (left.num * right.denom) == (left.denom * right.num);
  }

  friend bool operator<(const Rational& left, const Rational& right) {
    return (left.num * right.denom) < (left.denom * right.num);
  }

  bool operator!=(const Rational& right) const {
    return !(*this == right);
  }

  Rational& operator+=(const Rational& other) {
    if (denom == other.denom) {
      num += other.num;
      return *this;
    }

    num *= other.denom;
    num += (other.num * denom);
    denom *= other.denom;
    reduct();

    return *this;
  }

  Rational& operator-=(const Rational& other) {
    if (denom == other.denom) {
      num -= other.num;
      return *this;
    }

    num *= other.denom;
    num -= (other.num * denom);
    denom *= other.denom;
    reduct();

    return *this;
  }

  Rational& operator*=(const Rational& other) {
    num *= other.num;
    denom *= other.denom;
    reduct();

    return *this;
  }

  Rational& operator/=(const Rational& other) {
    if (*this == 0) {
      return *this;
    }
    if (other.num == 0) {
      throw std::exception();
    }

    num *= other.denom;
    denom *= other.num;
    if (denom.getSign() == num.getSign()) {
      num.getSign() = 1;
    } else {
      num.getSign() = -1;
    }
    denom.getSign() = 1;
    reduct();

    return *this;
  }

  Rational operator-() const {
    Rational res = *this;
    res.num.getSign() = (res.num.getSign() == 1) ? -1 : 1;

    return res;
  }

  explicit operator double() const {
    std::istringstream iss(asDecimal(kCastToDoublePrecision));

    double res;
    iss >> res;

    return res;
  }

  int& getSign() {
    return num.getSign();
  }

  const int& getSign() const {
    return num.getSign();
  }

  friend std::istream& operator>>(std::istream& in, Rational& real) {
    in >> real.num;
    real.denom = 1;

    return in;
  }

  friend std::ostream& operator<<(std::ostream& out, const Rational real) {
    out << real.toString();

    return out;
  }
};

Rational operator+(const Rational& left, const Rational& right) {
  Rational res = left;

  res += right;

  return res;
}

Rational operator-(const Rational& left, const Rational& right) {
  Rational res = left;

  res -= right;

  return res;
}

Rational operator*(const Rational& left, const Rational& right) {
  Rational res = left;

  res *= right;

  return res;
}

Rational operator/(const Rational& left, const Rational& right) {
  Rational res = left;

  res /= right;

  return res;
}

bool operator<=(const Rational& left, const Rational& right) {
  return (left < right || left == right);
}

bool operator>(const Rational& left, const Rational& right) {
  return !(left <= right);
}

bool operator>=(const Rational& left, const Rational& right) {
  return !(left < right);
}

Rational abs(Rational val) {
  val.getSign() = 1;

  return val;
}

template <unsigned long V>
struct Sqrt {
  template <int64_t L, int64_t R>
  struct SqrtHelper {
    static const int64_t M = (L + R + 1) / 2;
    static const int64_t N = R * R > V ? L : R;
    static const int64_t value = SqrtHelper<L + 1 >= R ? N : (M * M < V ? M : L),
                                            L + 1 >= R ? N : (M * M < V ? R : M)>::value;
  };
  template <int64_t M>
  struct SqrtHelper<M, M> {
    static const int64_t value = M;
  };
  static const int64_t value = SqrtHelper<0, V>::value;
};

template <int64_t N, int64_t K>
struct IsPrimeHelper {
  static const bool value = N % K != 0 && IsPrimeHelper<N, K - 1>::value;
};

template <int64_t N>
struct IsPrimeHelper<N, 2> {
  static const bool value = N % 2 != 0;
};

template <int64_t N>
struct IsPrimeHelper<N, 1> {
  static const bool value = N % 2 != 0;
};

template <unsigned long V>
struct IsPrime {
  static const bool value = IsPrimeHelper<V, Sqrt<V>::value>::value;
};

template <>
struct IsPrime<2> {
  static const bool value = true;
};

template <size_t N>
class Residue {
  int value = 0;

  static int ExtendedGCD(int a, int b, int& x, int& y) {
    if (a == 0) {
      x = 0;
      y = 1;
      return b;
    }
    int x1, y1;
    int d = ExtendedGCD(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
  }

  [[nodiscard]] int getInvertedValue() const {
    int u, v;
    ExtendedGCD(value, static_cast<int>(N), u, v);

    int res_v = normaliseValue(u);

    return normaliseValue(res_v);
  }

 public:
  Residue() = default;

  static int normaliseValue(const int _value) {
    return (_value % static_cast<int>(N) + static_cast<int>(N)) % static_cast<int>(N);
  }

  explicit Residue(const int other) {
    value = normaliseValue(other);
  }

  explicit operator int() const {
    return value;
  }

  int Value() const {
    return value;
  }

  Residue<N>& operator+=(const Residue<N>& other) {
    value += other.value;
    value = normaliseValue(value);

    return *this;
  }

  Residue<N>& operator-=(const Residue<N>& other) {
    value -= other.value;
    value = normaliseValue(value);

    return *this;
  }

  Residue<N>& operator*=(const Residue<N>& other) {
    value *= other.value;
    value = normaliseValue(value);

    return *this;
  }

  Residue<N>& operator/=(const Residue<N>& other) {
    static_assert(IsPrime<N>::value, "Division in residue field with non-prime base.");

    int inv_v = other.getInvertedValue();

    value *= inv_v;
    value = normaliseValue(value);

    return *this;
  }

  Residue<N> operator-() const {
    Residue<N> res = *this;
    res.value = normaliseValue(-res.value);

    return res;
  }

  friend std::ostream& operator<<(std::ostream& out, const Residue<N> rsd) {
    out << rsd.value;

    return out;
  }
};

template <size_t N>
Residue<N> operator+(const Residue<N>& left, const Residue<N>& right) {
  Residue<N> new_val = left;
  new_val += right;

  return new_val;
}

template <size_t N>
Residue<N> operator-(const Residue<N>& left, const Residue<N>& right) {
  Residue<N> new_val = left;
  new_val -= right;

  return new_val;
}

template <size_t N>
Residue<N> operator*(const Residue<N>& left, const Residue<N>& right) {
  Residue<N> new_val = left;
  new_val *= right;

  return new_val;
}

template <size_t N>
Residue<N> operator/(const Residue<N>& left, const Residue<N>& right) {
  Residue<N> new_val = left;
  new_val /= right;

  return new_val;
}

template <size_t N>
bool operator==(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() == right.Value();
}

template <size_t N>
bool operator<(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() < right.Value();
}

template <size_t N>
bool operator>(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() > right.Value();
}

template <size_t N>
bool operator<=(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() <= right.Value();
}

template <size_t N>
bool operator>=(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() >= right.Value();
}

template <size_t N>
bool operator!=(const Residue<N>& left, const Residue<N>& right) {
  return left.Value() != right.Value();
}

template <size_t N>
Residue<N> abs(const Residue<N> obj) {
  return obj;
}


template <size_t M, size_t N, typename Field=Rational>
class Matrix {
 public:
  static constexpr int IsCorrectHelper(const std::initializer_list<std::initializer_list<int>>& list) {
    if (list.size() != M) {
      return 1;
    }
    for (auto sublist : list) {
      if (sublist.size() != N) {
        return 1;
      }
    }

    return 0;
  }

  std::vector<std::vector<Field>> matrix_data = std::vector<std::vector<Field>>(M, std::vector<Field>(N, Field(0)));

  template <typename T>
  Matrix(std::initializer_list<std::initializer_list<T>> lst) {
    static_assert(IsCorrectHelper(lst) == 0, "Incorrect matrix size.");
    int i = 0, j = 0;
    for (auto& rw : lst) {
      for (auto& obj : rw) {
        matrix_data[i][j] = Field(obj);
        ++j;
      }
      ++i;
      j = 0;
    }
  }

  Matrix() = default;

  const std::vector<Field>& operator[](const size_t index) const {
    return matrix_data[index];
  }

  std::vector<Field>& operator[](const size_t index) {
    return matrix_data[index];
  }

  Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& other) {
    for (int i = 0; i < static_cast<int>(M); ++i) {
      for (int j = 0; j < static_cast<int>(N); ++j) {
        matrix_data[i][j] += other.matrix_data[i][j];
      }
    }

    return *this;
  }

  Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& other) {
    for (int i = 0; i < static_cast<int>(M); ++i) {
      for (int j = 0; j < static_cast<int>(N); ++j) {
        matrix_data[i][j] -= other.matrix_data[i][j];
      }
    }

    return *this;
  }

  Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& other) const {
    Matrix<M, N, Field> res = *this;
    res -= other;

    return res;
  }

  Matrix<M, N, Field>& operator*=(const Field other) {
    for (int i = 0; i < static_cast<int>(M); ++i) {
      for (int j = 0; j < static_cast<int>(N); ++j) {
        matrix_data[i][j] *= other;
      }
    }

    return *this;
  }

  Matrix<M, N, Field>& operator/=(const Field other) {
    *this *= (Field(1) / other);

    return *this;
  }

  Matrix<M, N, Field> operator*(const Field other) const {
    Matrix<M, N, Field> res = *this;
    res *= other;

    return res;
  }

  Matrix<M, N, Field> operator/(const Field other) const {
    Matrix<M, N, Field> res = *this;
    res /= other;

    return res;
  }

  template <size_t K>
  Matrix<M, N, Field>& operator*=(const Matrix<N, K, Field>& other) {
    static_assert(N == K, "Matrix dimensions must be equal");
    *this = (*this) * other;

    return *this;
  }

 private:
  static void swap(std::vector<std::vector<Field>> vec_c, int row1, int row2, int col) {
    for (int i = 0; i < col; ++i) {
      std::swap(vec_c[row1][i], vec_c[row2][i]);
    }
  }

 public:
  Field det() const {
    static_assert(M == N, "Determinant can be calculated only for square matrix.");

    std::vector<std::vector<Field>> vec_c = matrix_data;
    Field val = Field(1);

    for (int i = 0; i < static_cast<int>(M); ++i) {
      int k = i;
      for (int j = i + 1; j < static_cast<int>(M); ++j) {
        if (abs(vec_c[j][i]) > abs(vec_c[k][i])) {
          k = j;
        }
      }
      if (vec_c[k][i] == Field(0)) {
        val = Field(0);
        break;
      }

      std::swap(vec_c[i], vec_c[k]);
      if (i != k) {
        val = -val;
      }
      val *= vec_c[i][i];
      for (int j = i + 1; j < static_cast<int>(M); ++j) {
        vec_c[i][j] /= vec_c[i][i];
      }
      for (int j = 0; j < static_cast<int>(M); ++j) {
        if (j != i && vec_c[j][i] != Field(0)) {
          for (int l = i + 1; l < static_cast<int>(M); ++l) {
            vec_c[j][l] -= vec_c[i][l] * vec_c[j][i];
          }
        }
      }
    }

    return val;
  }

  Matrix<N, M, Field> transposed() const {
    Matrix<N, M, Field> matrixToTranspose;
    for (int i = 0; i < static_cast<int>(M); ++i) {
      for (int j = 0; j < static_cast<int>(N); ++j) {
        matrixToTranspose[j][i] = matrix_data[i][j];
      }
    }

    return matrixToTranspose;
  }

  [[nodiscard]] int rank() const {
    std::vector<std::vector<Field>> vec_c = matrix_data;
    int rank = static_cast<int>(N);

    for (int row = 0; row < rank; ++row) {
      if (vec_c[row][row] != Field(0)) {
        for (int col = 0; col < static_cast<int>(M); ++col) {
          if (col != row) {
            Field mult = vec_c[col][row] / vec_c[row][row];

            for (int i = 0; i < rank; ++i) {
              vec_c[col][i] -= mult * vec_c[row][i];
            }
          }
        }
      } else {
        bool reduce = true;

        for (int i = row + 1; i < static_cast<int>(M); ++i) {
          if (vec_c[i][row] != Field(0)) {
            swap(vec_c, row, i, rank);
            reduce = false;
            break;
          }
        }

        if (reduce) {
          --rank;

          for (int i = 0; i < static_cast<int>(M); ++i) {
            vec_c[i][row] = vec_c[i][rank];
          }
        }

        --row;
      }
    }

    return rank;
  }

  void invert() {
    static_assert(N == M, "Only square matrix can be inverted.");

    *this = this->inverted();
  }

  Matrix<M, M, Field> inverted() const {
    static_assert(N == M, "Only square matrix can be inverted.");
    int M_ = static_cast<int>(M);
    std::vector<std::vector<Field>> vec_c(M, std::vector<Field>(2 * M, Field(0)));
    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j < M_; ++j) {
        vec_c[i][j] = matrix_data[i][j];
      }
    }

    Field val;
    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j < M_; ++j) {
        if (i == j) {
          vec_c[i][j + M_] = Field(1);
        } else {
          vec_c[i][j + M_] = Field(0);
        }
      }
    }

    for (int i = 0; i < M_; ++i) {
      for (int j = 0; j < M_; ++j) {
        if (i != j) {
          val = vec_c[j][i] / vec_c[i][i];
          for (int k = 0; k < 2 * M_; ++k) {
            vec_c[j][k] -= val * vec_c[i][k];
          }
        }

      }
    }

    for (int i = 0; i < M_; ++i) {
      for (int j = M_; j < 2 * M_; ++j) {
        vec_c[i][j] /= vec_c[i][i];
      }
    }

    Matrix<M, M, Field> ans;
    for (int i = 0; i < M_; ++i) {
      for (int j = M_; j < 2 * M_; ++j) {
        ans[i][j - M_] = vec_c[i][j];
      }
    }

    return ans;
  }

  Field trace() const {
    static_assert(N == M, "Trace can be calculated only for square matrix.");

    Field val = Field(0);
    for (int i = 0; i < static_cast<int>(M); ++i) {
      val += matrix_data[i][i];
    }

    return val;
  }

  std::vector<Field> getRow(size_t row) {
    return matrix_data[row];
  }

  std::vector<Field> getColumn(size_t col) {
    std::vector<Field> columnVec(N);
    for (int i = 0; i < static_cast<int>(N); ++i) {
      columnVec[i] = matrix_data[col][i];
    }

    return columnVec;
  }

  friend std::ostream& operator<<(std::ostream& out, Matrix<M, N, Field> mtr) {
    for (int i = 0; i < static_cast<int>(M); ++i) {
      for (int j = 0; j < static_cast<int>(N); ++j) {
        out << mtr[i][j] << " ";
      }
      out << "\n";
    }

    return out;
  }
};

template <size_t N, size_t M, typename Field=Rational>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>* left, const Matrix<M, N, Field>& right) {
  Matrix<M, N, Field> res = left;
  res += right;

  return res;
}

template <size_t K, size_t N, size_t M, typename Field=Rational>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& left, const Matrix<N, K, Field>& right) {
  Matrix<M, K, Field> res;
  for (int i = 0; i < static_cast<int>(M); ++i) {
    for (int j = 0; j < static_cast<int>(K); ++j) {
      Field val = Field(0);
      for (int k = 0; k < static_cast<int>(N); ++k) {
        val += left[i][k] * right[k][j];
      }
      res[i][j] = val;
    }
  }

  return res;
}

template <size_t M, size_t N, typename Field=Rational>
bool operator==(const Matrix<M, N, Field>& left, const Matrix<M, N, Field>& right) {
  for (int i = 0; i < static_cast<int>(M); ++i) {
    for (int j = 0; j < static_cast<int>(N); ++j) {
      if (left[i][j] != right[i][j]) {
        return false;
      }
    }
  }
  return true;
}

template <size_t M, size_t N, typename Field=Rational>
bool operator!=(const Matrix<M, N, Field>& left, const Matrix<M, N, Field>& right) {
  return !(left == right);
}

template <size_t M, size_t N, typename Field=Rational>
Matrix<M, N, Field> operator*(const Field left, const Matrix<M, N, Field> right) {
  return right * left;
}

template <size_t N, typename Field=Rational>
using SquareMatrix = Matrix<N, N, Field>;

#endif //MATRIX__MATRIX_H_
