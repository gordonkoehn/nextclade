#pragma once


#include <fmt/format.h>


class Logger {
public:
  enum class Verbosity : int {
    silent = 0,
    error = 1,
    warn = 2,
    info = 3,
    debug = 4,
  };

  struct Options {
    Verbosity verbosity = Verbosity::warn;
  };


private:
  Options options;


public:
  Logger() = default;

  explicit Logger(const Options& loggerOptions) : options(loggerOptions) {}

  template<typename S, typename... Args>
  inline void debug(const S& format_str, Args&&... args) {
    if (options.verbosity >= Verbosity::debug) {
      fmt::print(format_str, args...);
    }
  }

  template<typename S, typename... Args>
  inline void info(const S& format_str, Args&&... args) {
    if (options.verbosity >= Verbosity::info) {
      fmt::print(format_str, args...);
    }
  }

  template<typename S, typename... Args>
  inline void warn(const S& format_str, Args&&... args) {
    if (options.verbosity >= Verbosity::warn) {
      fmt::print(format_str, args...);
    }
  }

  template<typename S, typename... Args>
  inline void error(const S& format_str, Args&&... args) {
    if (options.verbosity >= Verbosity::error) {
      fmt::print(stderr, format_str, args...);
    }
  }
};
