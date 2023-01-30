
#ifndef PERF_PH_UTIL_H
#define PERF_PH_UTIL_H

namespace tinspin {

template <typename Point>
class Box {
  public:
    explicit Box() = default;

    Box(const Box<Point>& orig) = default;

    Box(const Point& min, const Point& max) : min_{min}, max_{max} {}

    [[nodiscard]] const Point& min() const {
        return min_;
    }

    [[nodiscard]] const Point& max() const {
        return max_;
    }

    [[nodiscard]] Point& min() {
        return min_;
    }

    [[nodiscard]] Point& max() {
        return max_;
    }

    void min(const Point& new_min) {
        min_ = new_min;
    }

    void max(const Point& new_max) {
        max_ = new_max;
    }

    auto operator==(const Box<Point>& other) const -> bool {
        return min_ == other.min_ && max_ == other.max_;
    }

    auto operator!=(const Box<Point>& other) const -> bool {
        return !(*this == other);
    }

  private:
    Point min_;
    Point max_;
};
}  // namespace tinspin

#endif  // PERF_PH_UTIL_H
