use std::cmp::Reverse;
use std::collections::BinaryHeap;

pub const INF: u64 = u64::MAX / 4;

#[derive(Clone, Debug)]
pub struct Graph {
    num_nodes: usize,
    edges: Vec<Vec<(usize, u64)>>,
}

impl Graph {
    pub fn new(num_nodes: usize) -> Self {
        Self { num_nodes, edges: vec![Vec::new(); num_nodes] }
    }

    pub fn num_nodes(&self) -> usize {
        self.num_nodes
    }

    pub fn add_edge(&mut self, from: usize, to: usize, weight: u64) {
        assert!(from < self.num_nodes && to < self.num_nodes);
        self.edges[from].push((to, weight));
    }

    pub fn add_undirected_edge(&mut self, u: usize, v: usize, weight: u64) {
        self.add_edge(u, v, weight);
        self.add_edge(v, u, weight);
    }

    pub fn neighbors(&self, u: usize) -> &[(usize, u64)] {
        &self.edges[u]
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Distances {
    pub distances: Vec<u64>,
    pub predecessors: Vec<usize>,
}

impl Distances {
    pub fn path_to(&self, target: usize) -> Option<Vec<usize>> {
        if target >= self.distances.len() || self.distances[target] == INF {
            return None;
        }
        let mut path_len = 1usize;
        let mut cur = target;
        while self.predecessors[cur] != cur {
            if self.predecessors[cur] == usize::MAX {
                return None;
            }
            path_len += 1;
            cur = self.predecessors[cur];
        }
        let mut path = vec![0usize; path_len];
        let mut i = path_len - 1;
        let mut cur = target;
        loop {
            path[i] = cur;
            if self.predecessors[cur] == cur { break; }
            cur = self.predecessors[cur];
            if i == 0 { break; }
            i -= 1;
        }
        Some(path)
    }
}

pub fn bounded_multi_source_shortest_paths(graph: &Graph, sources: &[usize], bound: u64) -> Distances {
    let n = graph.num_nodes();
    for &s in sources { debug_assert!(s < n); }
    let mut dist = vec![INF; n];
    let mut pred = vec![usize::MAX; n];
    let mut heap: BinaryHeap<(Reverse<u64>, usize)> = BinaryHeap::new();
    for &s in sources {
        if dist[s] != 0 {
            dist[s] = 0;
            pred[s] = s;
            heap.push((Reverse(0), s));
        }
    }
    while let Some((Reverse(d_u), u)) = heap.pop() {
        if d_u > bound { break; }
        if d_u != dist[u] { continue; }
        for &(v, w) in graph.neighbors(u) {
            if w <= bound.saturating_sub(d_u) {
                let nd = d_u + w;
                if nd < dist[v] {
                    dist[v] = nd;
                    pred[v] = u;
                    heap.push((Reverse(nd), v));
                }
            }
        }
    }
    Distances { distances: dist, predecessors: pred }
}

pub fn unbounded_multi_source_shortest_paths(graph: &Graph, sources: &[usize]) -> Distances {
    let n = graph.num_nodes();
    for &s in sources { debug_assert!(s < n); }
    let mut dist = vec![INF; n];
    let mut pred = vec![usize::MAX; n];
    let mut heap: BinaryHeap<(Reverse<u64>, usize)> = BinaryHeap::new();
    for &s in sources {
        if dist[s] != 0 {
            dist[s] = 0;
            pred[s] = s;
            heap.push((Reverse(0), s));
        }
    }
    while let Some((Reverse(d_u), u)) = heap.pop() {
        if d_u != dist[u] { continue; }
        for &(v, w) in graph.neighbors(u) {
            let nd = d_u.saturating_add(w);
            if nd < dist[v] {
                dist[v] = nd;
                pred[v] = u;
                heap.push((Reverse(nd), v));
            }
        }
    }
    Distances { distances: dist, predecessors: pred }
}

#[cfg(test)]
mod tests {
    use super::*;

    

    #[test]
    fn line_graph_bounded() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1);
        g.add_edge(1, 2, 1);
        g.add_edge(2, 3, 1);
        let d = bounded_multi_source_shortest_paths(&g, &[0], 2);
        assert_eq!(d.distances[0], 0);
        assert_eq!(d.distances[1], 1);
        assert_eq!(d.distances[2], 2);
        assert_eq!(d.distances[3], INF);
        let p = d.path_to(2).unwrap();
        assert_eq!(p, vec![0,1,2]);
        assert!(d.path_to(3).is_none());
    }

    #[test]
    fn multi_source_star_bound_1() {
        let mut g = Graph::new(5);
        for v in 1..5 { g.add_undirected_edge(0, v, 1); }
        let d = bounded_multi_source_shortest_paths(&g, &[1,3], 1);
        assert_eq!(d.distances[1], 0);
        assert_eq!(d.distances[3], 0);
        assert_eq!(d.distances[0], 1);
        assert_eq!(d.distances[2], INF);
        assert_eq!(d.distances[4], INF);
    }

    #[test]
    fn multi_source_star_bound_2() {
        let mut g = Graph::new(5);
        for v in 1..5 { g.add_undirected_edge(0, v, 1); }
        let d = bounded_multi_source_shortest_paths(&g, &[1,3], 2);
        assert_eq!(d.distances[1], 0);
        assert_eq!(d.distances[3], 0);
        assert_eq!(d.distances[0], 1);
        assert_eq!(d.distances[2], 2);
        assert_eq!(d.distances[4], 2);
    }

    #[test]
    fn grid_3x3_bounded() {
        let n = 3;
        let idx = |r: usize, c: usize| -> usize { r * n + c };
        let mut g = Graph::new(n * n);
        for r in 0..n {
            for c in 0..n {
                let u = idx(r, c);
                if r + 1 < n { g.add_undirected_edge(u, idx(r + 1, c), 1); }
                if c + 1 < n { g.add_undirected_edge(u, idx(r, c + 1), 1); }
            }
        }
        let d = bounded_multi_source_shortest_paths(&g, &[idx(0,0), idx(2,2)], 2);
        assert_eq!(d.distances[idx(1,1)], 2);
        assert!(d.distances.iter().all(|&x| x <= 2 || x == INF));
    }

    #[test]
    fn zero_weight_edges() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 0);
        g.add_edge(1, 2, 0);
        g.add_edge(2, 3, 0);
        let d = bounded_multi_source_shortest_paths(&g, &[0], 0);
        assert_eq!(d.distances[0], 0);
        assert_eq!(d.distances[1], 0);
        assert_eq!(d.distances[2], 0);
        assert_eq!(d.distances[3], 0);
        let p = d.path_to(3).unwrap();
        assert_eq!(p, vec![0,1,2,3]);
    }

    #[test]
    fn disconnected_graph() {
        let mut g = Graph::new(5);
        g.add_edge(0, 1, 2);
        g.add_edge(1, 2, 2);
        g.add_edge(3, 4, 1);
        let d = bounded_multi_source_shortest_paths(&g, &[0], 10);
        assert_eq!(d.distances[3], INF);
        assert_eq!(d.distances[4], INF);
    }

    #[test]
    fn bound_zero_only_sources() {
        let mut g = Graph::new(4);
        g.add_edge(0, 1, 1);
        g.add_edge(1, 2, 1);
        g.add_edge(2, 3, 1);
        let d = bounded_multi_source_shortest_paths(&g, &[1,3], 0);
        assert_eq!(d.distances[0], INF);
        assert_eq!(d.distances[1], 0);
        assert_eq!(d.distances[2], INF);
        assert_eq!(d.distances[3], 0);
    }

    #[test]
    fn bounded_matches_unbounded_with_mask() {
        let mut g = Graph::new(6);
        g.add_undirected_edge(0, 1, 1);
        g.add_undirected_edge(1, 2, 3);
        g.add_undirected_edge(2, 3, 1);
        g.add_undirected_edge(3, 4, 4);
        g.add_undirected_edge(1, 5, 2);
        let ub = unbounded_multi_source_shortest_paths(&g, &[0,4]);
        let b = bounded_multi_source_shortest_paths(&g, &[0,4], 3);
        for i in 0..g.num_nodes() {
            if ub.distances[i] <= 3 {
                assert_eq!(b.distances[i], ub.distances[i]);
            } else {
                assert_eq!(b.distances[i], INF);
            }
        }
    }

    #[test]
    fn duplicate_sources_handled() {
        let mut g = Graph::new(3);
        g.add_undirected_edge(0, 1, 1);
        g.add_undirected_edge(1, 2, 1);
        let d = bounded_multi_source_shortest_paths(&g, &[0,0], 2);
        assert_eq!(d.distances, vec![0,1,2]);
    }

    #[test]
    fn large_weights_no_overflow() {
        let mut g = Graph::new(3);
        g.add_edge(0, 1, u64::MAX/8);
        g.add_edge(1, 2, u64::MAX/8);
        let d = bounded_multi_source_shortest_paths(&g, &[0], u64::MAX/4);
        assert!(d.distances[1] < INF);
        assert!(d.distances[2] < INF);
    }
}
