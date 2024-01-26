#![allow(dead_code)]
use identity_hash::{IdentityHashable, IntMap, IntSet};
use std::collections::HashMap;

pub type NodeID = usize;

pub struct Node<T: IdentityHashable, S> {
    pub children: Vec<NodeID>,
    content: Vec<T>,
    pub parent: Option<NodeID>,
    pub tag: S,
    pub metadata: HashMap<String, String>,
}

impl<T: IdentityHashable, S> Node<T, S> {
    pub fn content_is_empty(&self) -> bool {
        self.content.is_empty()
    }
    pub fn content_len(&self) -> usize {
        self.content.len()
    }
    pub fn content_iter(&self) -> impl Iterator<Item = &T> {
        self.content.iter()
    }
    pub fn content_slice(&self) -> &[T] {
        &self.content
    }
    pub fn content_clear(&mut self) {
        self.content.clear()
    }
}

pub struct PTree<T: Clone + IdentityHashable, S> {
    current_id: usize,
    nodes: HashMap<NodeID, Node<T, S>>,
    descendants_cache: IntMap<NodeID, IntSet<NodeID>>,
    leaves_cache: IntMap<NodeID, IntSet<T>>,
    leaves_cache_enabled: bool,
}

impl<T: Clone + IdentityHashable, S> std::ops::Index<usize> for PTree<T, S> {
    type Output = Node<T, S>;
    fn index(&self, i: usize) -> &Self::Output {
        &self.nodes[&i]
    }
}
impl<T: Clone + IdentityHashable, S> std::ops::IndexMut<usize> for PTree<T, S> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        self.nodes.get_mut(&i).unwrap()
    }
}

impl<T: Clone + Eq + IdentityHashable + std::hash::Hash + std::fmt::Debug, S> PTree<T, S> {
    pub fn new() -> Self {
        PTree {
            current_id: 0,
            nodes: HashMap::new(),
            descendants_cache: Default::default(),
            leaves_cache: Default::default(),
            leaves_cache_enabled: false,
        }
    }

    fn add_to_cache(&mut self, n: NodeID, new: impl IntoIterator<Item = NodeID> + Clone) {
        if let Some(descs) = self.descendants_cache.get_mut(&n) {
            for n in new.clone() {
                descs.insert(n);
            }
        }
        if let Some(parent) = self[n].parent {
            self.add_to_cache(parent, new.clone());
        }
    }

    fn rm_from_cache(&mut self, n: NodeID, old: impl IntoIterator<Item = NodeID> + Clone) {
        if let Some(descs) = self.descendants_cache.get_mut(&n) {
            for n in old.clone() {
                descs.remove(&n);
            }
        }
        if let Some(parent) = self[n].parent {
            self.rm_from_cache(parent, old.clone());
        }
    }

    fn add_to_cache_and_children(&mut self, n: NodeID) {
        if let Some(parent) = self[n].parent {
            self.add_to_cache(parent, [n]);
            self.add_to_cache(
                parent,
                self.descendants_cache[&n].iter().cloned().collect::<Vec<_>>(),
            );
        }
    }

    fn rm_from_cache_and_children(&mut self, n: NodeID) {
        if let Some(parent) = self[n].parent {
            self.rm_from_cache(parent, [n]);
            self.rm_from_cache(
                parent,
                self.descendants_cache[&n].iter().cloned().collect::<Vec<_>>(),
            );
        }
    }

    fn cache_leaves(&mut self, n: NodeID, xs: &[T]) {
        for x in xs {
            self.cache_leave(n, x.clone());
        }
    }

    fn cache_leave(&mut self, n: NodeID, x: T) {
        if !self.leaves_cache_enabled {
            return;
        }

        self.leaves_cache.get_mut(&n).unwrap().insert(x.clone());

        if let Some(parent) = self[n].parent {
            self.cache_leave(parent, x);
        }
    }

    // fn uncache_leave(&mut self, n: NodeID, x: T) {
    //     if let Some(parent) = self[n].parent {
    //         self.leaves_cache.get_mut(&parent).unwrap().get_mut(&x).map(|k| {
    //             if *k > 0 {
    //                 *k -= 1
    //             }
    //         });
    //         self.uncache_leave(parent, x);
    //     }
    // }

    // fn uncache_leaves(&mut self, n: NodeID, xs: &[T]) {
    //     for x in xs {
    //         self.uncache_leave(n, x.clone());
    //     }
    // }

    pub fn add_leave(&mut self, n: NodeID, x: T) {
        self[n].content.push(x.clone());
        self.cache_leave(n, x);
    }

    pub fn clear_leaves(&mut self, n: NodeID) {
        self[n].content.clear();
        self.refresh_cache_leaves(n);
    }

    fn rec_refresh_cache_leaves(&mut self, n: NodeID) {
        if let Some(parent) = self[n].parent {
            *self.leaves_cache.get_mut(&parent).unwrap() = self[parent].children.iter().fold(
                self[parent].content.iter().cloned().collect::<IntSet<_>>(),
                |mut ax, c| {
                    for cc in self.leaves_cache[c].iter().cloned() {
                        ax.insert(cc);
                    }
                    ax
                },
            );
            self.rec_refresh_cache_leaves(parent)
        }
    }

    fn refresh_cache_leaves(&mut self, n: NodeID) {
        if !self.leaves_cache_enabled {
            return;
        }
        *self.leaves_cache.get_mut(&n).unwrap() = self.full_descendant_leaves(n);
        self.rec_refresh_cache_leaves(n);
    }

    fn pure_refresh_cache_leaves(&mut self, n: NodeID) {
        if !self.leaves_cache_enabled {
            return;
        }
        self.rec_refresh_cache_leaves(n);
    }

    pub fn disable_leave_cache(&mut self) {
        self.leaves_cache_enabled = false;
        self.leaves_cache.clear();
    }

    pub fn enable_leave_cache(&mut self) {
        self.leaves_cache_enabled = true;
        let nodes = self.nodes().cloned().collect::<Vec<_>>();

        for n in nodes.into_iter() {
            self.leaves_cache.insert(n, self.full_descendant_leaves(n));
        }
        // for l in self.nodes().filter(|n| self[**n].children.is_empty()) {
        //     self.leaves_cache.insert(*l, self[*l].content.iter().cloned().collect());
        // }

        // for l in self.nodes().filter(|n| self[**n].children.is_empty()) {
        //     self.pure_refresh_cache_leaves(*l);
        // }
    }

    pub fn add_node(&mut self, content: &[T], tag: S, parent: Option<NodeID>) -> NodeID {
        self.current_id = self.current_id.checked_add(1).expect("Tree is too big");
        let id = self.current_id;
        assert!(!self.nodes.contains_key(&id), "{} already exists", id);
        assert!(parent.is_none() || self.nodes.contains_key(&parent.unwrap()));
        self.nodes.insert(
            id,
            Node {
                children: Vec::with_capacity(2),
                content: content.to_vec(),
                parent,
                tag,
                metadata: Default::default(),
            },
        );
        self.descendants_cache.insert(id, Default::default());
        self.leaves_cache.insert(id, Default::default());
        if let Some(parent) = parent {
            self[parent].children.push(id);
            self.add_to_cache(parent, [id]);
        }
        self.cache_leaves(id, content);
        id
    }

    pub fn nodes(&self) -> impl Iterator<Item = &NodeID> {
        self.nodes.keys()
    }

    pub fn plug(&mut self, target: NodeID, n: NodeID) {
        assert!(self.nodes[&n].parent.is_none());
        assert!(!self.nodes[&target].children.contains(&n));
        self.nodes.get_mut(&n).unwrap().parent = Some(target);
        self.nodes.get_mut(&target).unwrap().children.push(n);
        self.add_to_cache_and_children(n);
        self.pure_refresh_cache_leaves(n);
    }

    pub fn unplug(&mut self, n: NodeID) {
        let parent = self.nodes[&n].parent;
        assert!(parent.is_none() || self.nodes[&parent.unwrap()].children.contains(&n));

        if let Some(parent) = parent {
            self.nodes.get_mut(&parent).unwrap().children.retain(|nn| *nn != n);
            self.rm_from_cache_and_children(n);
            self.pure_refresh_cache_leaves(n);
        }
        self.nodes.get_mut(&n).unwrap().parent = None;
    }

    pub fn unplug_many(&mut self, parent: NodeID, ns: &[NodeID]) {
        for &n in ns {
            assert!(self.nodes[&parent].children.contains(&n));
            self.unplug(n);
        }
    }

    pub fn delete_node(&mut self, n: NodeID) {
        assert!(self.nodes.contains_key(&n));
        self.rm_from_cache_and_children(n);
        for c in self[n].children.clone().into_iter() {
            self.delete_node(c);
        }
        self.unplug(n);
        self.nodes.remove(&n);
    }

    pub fn delete_nodes(&mut self, ns: &[NodeID]) {
        ns.iter().for_each(|&n| self.delete_node(n));
    }

    pub fn move_node(&mut self, n: NodeID, dest: NodeID) {
        self.unplug(n);
        self.plug(dest, n);
    }

    fn rec_descendants(&self, i: NodeID, ax: &mut IntSet<NodeID>) {
        for j in &self[i].children {
            ax.insert(*j);
            self.rec_descendants(*j, ax)
        }
    }

    pub fn full_descendants(&self, n: NodeID) -> IntSet<NodeID> {
        let mut r = Default::default();
        self.rec_descendants(n, &mut r);
        r
    }

    pub fn descendants(&self, n: NodeID) -> &IntSet<NodeID> {
        &self.descendants_cache[&n]
    }

    pub fn topo_depth(&self, n: NodeID) -> usize {
        let mut me = n;
        let mut depth = 0;
        while let Some(parent) = self.nodes[&me].parent {
            depth += 1;
            me = parent;
        }
        depth
    }

    fn format_leaf_newick(
        &self,
        i: NodeID,
        f_leaf: &dyn Fn(&T) -> String,
        f_tag: &dyn Fn(&S) -> String,
    ) -> String {
        let mut r = String::new();

        let leaves = self.nodes[&i].content.iter().map(f_leaf).collect::<Vec<String>>().join(",");
        let children = self.nodes[&i]
            .children
            .iter()
            .map(|&c| self.format_leaf_newick(c, f_leaf, f_tag))
            .filter(|s| !s.is_empty())
            .collect::<Vec<String>>()
            .join(",");

        r.push('(');
        r.push_str(&leaves);
        if !leaves.is_empty() && !children.is_empty() {
            r.push(',');
        }
        r.push_str(&children);
        r.push(')');
        let metadata = self.nodes[&i]
            .metadata
            .iter()
            .map(|(k, v)| format!("{}={}", k, v))
            .collect::<Vec<String>>()
            .join(":");
        r.push_str(&format!(
            "[&&NHX:ID={i}:S={}{}]",
            f_tag(&self.nodes[&i].tag),
            if metadata.is_empty() { "".to_owned() } else { format!(":{}", metadata) }
        ));

        r
    }

    pub fn cardinal(&self, n: usize) -> usize {
        self[n].children.len() + self[n].content.len()
    }

    fn rec_descendant_leaves(&self, i: NodeID, ax: &mut IntSet<T>) {
        for l in &self.nodes[&i].content {
            ax.insert(l.clone());
        }
        // ax.extend_from_slice(&self.nodes[&i].content);
        for &j in &self.nodes[&i].children {
            self.rec_descendant_leaves(j, ax)
        }
    }
    pub fn full_descendant_leaves(&self, n: NodeID) -> IntSet<T> {
        let mut r = Default::default();
        self.rec_descendant_leaves(n, &mut r);
        r
    }

    pub fn descendant_leaves(&self, n: NodeID) -> &IntSet<T> {
        &self.leaves_cache[&n]
    }

    pub fn find_node<F>(&self, f: F) -> Option<NodeID>
    where
        F: Fn(&Node<T, S>) -> bool,
    {
        self.nodes.iter().find(|(_i, n)| f(n)).map(|(i, _)| *i)
    }

    pub fn to_newick(&self, f_leaf: &dyn Fn(&T) -> String, f_tag: &dyn Fn(&S) -> String) -> String {
        let mut r = String::new();

        for k in self.nodes.keys().filter(|k| self.nodes[k].parent.is_none()) {
            let mut k = *k;
            while self[k].children.len() == 1 && self[k].content.is_empty() {
                k = self[k].children[0];
            }
            r.push_str(&self.format_leaf_newick(k, f_leaf, f_tag));
            r.push_str(";\n");
        }

        r
    }

    pub fn to_newick_from(
        &self,
        subroot: usize,
        f_leaf: &dyn Fn(&T) -> String,
        f_tag: &dyn Fn(&S) -> String,
    ) -> String {
        let mut r = String::new();

        r.push_str(&self.format_leaf_newick(subroot, f_leaf, f_tag));
        r.push_str(";\n");

        r
    }
}

impl<T: Clone + Eq + std::fmt::Debug + IdentityHashable, S: std::fmt::Debug> PTree<T, S> {
    fn rec_disp(&self, i: usize, depth: usize) {
        println!("{}{}#{:?} > CHILDREN: {:?}", " ".repeat(depth), i, self[i].tag, self[i].children);
        println!("{}    |- CONTENT: {:?}", " ".repeat(depth), self[i].content);

        for j in self[i].children.iter() {
            self.rec_disp(*j, depth + 2);
        }
    }
    pub fn disp(&self, root: usize) {
        self.rec_disp(root, 0);
    }
}
