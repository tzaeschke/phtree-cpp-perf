
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

/**
 * Very simple but fast stack.
 *
 * WARNING: This stack violates expected behavior of normal stacks:
 * - It may construct more objects than required
 * - It may overwrite objects using the assignment operator
 * - Destructors may not be called (i.e. the object may be overwritten by assignment instead)
 * - ...
 */
template<typename T>
class SimpleStack {
  public:
    SimpleStack(size_t size = 10) : stack(size) {}

    bool empty() noexcept {
        return size == 0u;
    }

    T& push() noexcept {
        if (size == stack.size()) {
            stack.emplace_back();
        }
        return stack[size++];
    }

    template<typename... Args>
    T& push(Args&&... args) noexcept {
        if (size == stack.size()) {
            stack.emplace_back(std::forward<Args>(args)...);
        } else {
            // TODO this is obviously bad, we also never call the destructor.
            stack[size] = T{std::forward<Args>(args)...};
        }
        return stack[size++];
    }

    template<typename T2 = T>
    T& push(T&& t) noexcept {
        if (size == stack.size()) {
            stack.emplace_back(std::forward<T2>(t));
        } else {
            // TODO this is obviously bad, we also never call the destructor.
            stack[size] = T{std::forward<T2>(t)};
        }
        return stack[size++];
    }

    T& top() noexcept {
        return stack[size - 1];
    }

    T& pop() noexcept  {
        return stack[--size];
    }

    void clear() noexcept  {
        size = 0;
    }

  private:
    std::vector<T> stack{};
    size_t size = 0;
};


}  // namespace tinspin

#endif  // PERF_PH_UTIL_H
