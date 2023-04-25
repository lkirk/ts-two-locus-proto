#include <stdio.h>

int main() {
  int n = 3;
  int subloop_start = 0;
  for (int i = 0; i < n; i++) {
    for (int j = subloop_start; j < n; j++) {
      printf("%d\t%d\n", i, j);
    }
    subloop_start++;
  }
}
