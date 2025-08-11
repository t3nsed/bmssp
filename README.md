# Bounded Multi-Source Shortest Paths

Just a small fun excersise to implement the new bounded multi-source shortest paths problem in Rust.

We're looking at $O(m \log^{(2/3)} n)$ time. For the limited tests I've run, it's still like 10% slower than the current reference rust implementation of Dijkstra's SSSP, but I'm sure it's a matter of a few rust-specific compute and/or caching tweaks.

https://x.com/dorsa_rohani/status/1954573594853244964

https://arxiv.org/pdf/2504.17033