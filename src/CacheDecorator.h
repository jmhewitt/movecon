#ifndef MOVECON_CACHE_DECORATOR_H
#define MOVECON_CACHE_DECORATOR_H



#include <functional>
#include <map>

/**
 * source: https://martin-ueding.de/posts/c-cache-decorator/
*/
template <typename R, typename... A>
class CacheDecorator {
  public:
    CacheDecorator(std::function<R(A...)> f) : f_(f) {}

    R operator()(A... a) {
        std::tuple<A...> key(a...);
        auto search = map_.find(key);
        if (search != map_.end()) {
            return search->second;
        }

        auto result = f_(a...);
        map_[key] = result;
        return result;
    }

  private:
    std::function<R(A...)> f_;
    std::map<std::tuple<A...>, R> map_;
};

#endif
