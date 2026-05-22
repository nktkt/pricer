# 金融計算サンプル（C++ / LLVM）のビルド用 Makefile
CXX      := clang++
CXXFLAGS := -std=c++17 -O2
LLVM     := /opt/homebrew/opt/llvm/bin/llvm-config

# LLVM JIT 用のフラグ（-fexceptions で llvm-config の -fno-exceptions を上書き）
JITFLAGS := $(shell $(LLVM) --cxxflags --ldflags --libs core orcjit native --system-libs) -fexceptions

# 普通のサンプル（LLVM 不要）
PLAIN := option_pricing convergence parallel_mc greeks barrier_option

.PHONY: all run clean
all: $(PLAIN) jit_payoff

$(PLAIN): %: %.cpp bs_common.hpp
	$(CXX) $(CXXFLAGS) -pthread $< -o $@

jit_payoff: jit_payoff.cpp bs_common.hpp
	$(CXX) $(CXXFLAGS) $< $(JITFLAGS) -o $@

# 全サンプルを順に実行
run: all
	@for p in $(PLAIN); do echo "===== $$p ====="; ./$$p; echo; done
	@echo "===== jit_payoff ====="; ./jit_payoff

clean:
	rm -f $(PLAIN) jit_payoff
