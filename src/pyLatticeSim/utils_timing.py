import time
from collections import defaultdict
import threading
from colorama import Fore, Style


class Timing:
    """
    A class to time function execution and track call hierarchy.
    """
    def __init__(self):
        self.timings = defaultdict(list)
        self.call_stack = []  # to track call hierarchy
        self.call_graph = defaultdict(lambda: defaultdict(float))
        self.call_counts = defaultdict(int)
        self.local = threading.local()
        self.start_time = time.perf_counter()


    def timeit(self, func):
        """
        Decorator to time the execution of a function and track its calls.
        """
        def wrapper(*args, **kwargs):
            parent = self.call_stack[-1] if self.call_stack else None
            self.call_stack.append(func.__name__)
            start = time.perf_counter()
            result = func(*args, **kwargs)
            end = time.perf_counter()
            elapsed = end - start
            self.timings[func.__name__].append(elapsed)
            self.call_counts[func.__name__] += 1
            if parent:
                self.call_graph[parent][func.__name__] += elapsed
            self.call_stack.pop()
            return result
        return wrapper

    def summary(self):
        """
        Print a summary of the timings collected.
        """
        avg = lambda times: sum(times) / len(times) if times else 0
        print(Fore.GREEN + f"{'Function':<30} {'Calls':<10} {'Total (s)':<12} {'Avg (s)':<12} {'Max (s)':<12}")
        print("-" * 80 + Style.RESET_ALL)
        sorted_funcs = sorted(self.timings.items(), key=lambda x: sum(x[1]), reverse=True)
        for name, times in sorted_funcs:
            total = sum(times)
            count = len(times)
            print(f"{name:<50} {count:<10} {total:<12.6f} {avg(times):<12.6f} {max(times):<12.6f}")
            # If the function has callees
            if name in self.call_graph:
                for subname, subtime in self.call_graph[name].items():
                    print(f"  └─ {subname:<26} {self.call_counts[subname]:<10} {subtime:<12.6f}")

        total_elapsed = time.perf_counter() - self.start_time
        print(Fore.LIGHTYELLOW_EX, f"\n🔧 Total script runtime: {total_elapsed:.4f} s" + Style.RESET_ALL)