use std::collections::HashMap;

pub type NodeID = usize;
pub struct Node<T, S> {
    pub children: Vec<NodeID>,
    pub content: Vec<T>,
    pub parent: Option<NodeID>,
    pub tag: S,
}

pub struct PTree<T: Clone, S> {
    current_id: usize,
    nodes: HashMap<NodeID, Node<T, S>>,
    descendants_cache: HashMap<NodeID, Vec<NodeID>>,
}

impl<T: Clone, S> std::ops::Index<usize> for PTree<T, S> {
    type Output = Node<T, S>;
    fn index(&self, i: usize) -> &Self::Output {
        &self.nodes[&i]
    }
}
impl<T: Clone, S> std::ops::IndexMut<usize> for PTree<T, S> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        self.nodes.get_mut(&i).unwrap()
    }
}

impl<T: Clone, S> PTree<T, S> {
    pub fn new() -> Self {
        PTree {
            current_id: 0,
            nodes: HashMap::new(),
            descendants_cache: HashMap::new(),
        }
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
            },
        );
        if let Some(parent) = parent {
            self[parent].children.push(id);
        }
        id
    }

    pub fn nodes(&self) -> impl Iterator<Item = &NodeID> {
        self.nodes.keys()
    }

    pub fn plug(&mut self, target: NodeID, n: NodeID) {
        assert!(self.nodes[&n].parent.is_none());
        assert!(!self.nodes[&target].children.contains(&n));
        self.nodes.get_mut(&n).unwrap().parent = Some(target);
        self.nodes.get_mut(&target).unwrap().children.push(n)
    }

    pub fn unplug(&mut self, n: NodeID) {
        let parent = self.nodes[&n].parent;
        assert!(parent.is_none() || self.nodes[&parent.unwrap()].children.contains(&n));

        self.nodes.get_mut(&n).unwrap().parent = None;
        if let Some(parent) = parent {
            self.nodes
                .get_mut(&parent)
                .unwrap()
                .children
                .retain(|nn| *nn != n);
        }
    }

    pub fn unplug_many(&mut self, parent: NodeID, ns: &[NodeID]) {
        for &n in ns {
            assert!(self.nodes[&parent].children.contains(&n));
            self.unplug(n);
        }
    }

    pub fn delete_node(&mut self, n: NodeID) {
        assert!(self.nodes.contains_key(&n));
        for c in self[n].children.clone().into_iter() {
            self.delete_node(c);
        }
        self.unplug(n);
        self.nodes.remove(&n);
    }

    pub fn delete_nodes(&mut self, ns: &[NodeID]) {
        ns.iter().for_each(|&n| self.delete_node(n));
    }

    pub fn merge_nodes(
        &mut self,
        merger: NodeID,
        merged: NodeID,
        f: &dyn Fn(&mut Vec<T>, &[T]) -> (),
    ) {
        assert!(self.nodes.contains_key(&merger));
        assert!(self.nodes.contains_key(&merged));

        self.nodes
            .values_mut()
            .filter(|v| v.parent.is_some() && v.parent.unwrap() == merged)
            .for_each(|v| v.parent = Some(merger));

        let merged_children = self.nodes[&merged].children.to_vec();
        self.nodes
            .get_mut(&merger)
            .unwrap()
            .children
            .extend_from_slice(&merged_children);

        let merged_content = &self.nodes[&merged].content.clone();
        f(
            &mut self.nodes.get_mut(&merger).unwrap().content,
            &merged_content,
        );
        self.delete_node(merged);
    }

    pub fn move_node(&mut self, n: NodeID, dest: NodeID) {
        self.unplug(n);
        self.plug(dest, n);
    }

    fn rec_descendants(&self, i: NodeID, ax: &mut Vec<NodeID>) {
        for j in &self[i].children {
            ax.push(*j);
            self.rec_descendants(*j, ax)
        }
    }

    pub fn descendants(&self, n: NodeID) -> Vec<NodeID> {
        let mut r = Vec::with_capacity(self.nodes.len() / 2);
        self.rec_descendants(n, &mut r);
        r
    }

    pub fn cache_descendants(&mut self, from: NodeID) {
        let mut me = from;
        let todo = self.descendants(me);
        self.descendants_cache.insert(me, todo.to_owned());
        for n in todo {
            self.descendants_cache.insert(n, self.descendants(n));
        }
        while let Some(parent) = self[me].parent {
            self.descendants_cache
                .insert(parent, self.descendants(parent));
            me = parent;
        }
    }

    pub fn cached_descendants(&self, n: NodeID) -> Option<&Vec<NodeID>> {
        self.descendants_cache.get(&n)
    }

    fn rec_descendant_leaves(&self, i: NodeID, ax: &mut Vec<T>) {
        ax.extend_from_slice(&self.nodes[&i].content);
        for &j in &self.nodes[&i].children {
            self.rec_descendant_leaves(j, ax)
        }
    }
    pub fn descendant_leaves(&self, n: NodeID) -> Vec<T> {
        let mut r = Vec::with_capacity(self.nodes.len() / 2);
        self.rec_descendant_leaves(n, &mut r);
        r
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

        let ms = self.nodes[&i]
            .content
            .iter()
            .map(|c| f_leaf(c))
            .collect::<Vec<String>>()
            .join(",");
        let cs = self.nodes[&i]
            .children
            .iter()
            .map(|&c| self.format_leaf_newick(c, f_leaf, f_tag))
            .filter(|s| !s.is_empty())
            .collect::<Vec<String>>()
            .join(",");

        r.push_str("(");
        r.push_str(&ms);
        if !ms.is_empty() && !cs.is_empty() {
            r.push_str(",");
        }
        r.push_str(&cs);
        r.push_str(")");

        if !self.nodes[&i].children.is_empty() || self.nodes[&i].content.len() > 1 {
            r.push_str(&format!("{}-{}", f_tag(&self.nodes[&i].tag), i));
        }
        r.push_str(&format!("[&&NHX:S={}]", f_tag(&self.nodes[&i].tag)));

        r
    }
    pub fn to_newick(&self, f_leaf: &dyn Fn(&T) -> String, f_tag: &dyn Fn(&S) -> String) -> String {
        let mut r = String::new();

        for k in self
            .nodes
            .keys()
            .filter(|k| self.nodes[&k].parent.is_none())
        {
            let mut k = *k;
            while self[k].children.len() == 1 {
                k = self[k].children[0];
            }
            r.push_str(&self.format_leaf_newick(k, f_leaf, f_tag));
            r.push_str(";\n");
        }

        r
    }
    pub fn to_newick_from(
        &self,
        i: usize,
        f_leaf: &dyn Fn(&T) -> String,
        f_tag: &dyn Fn(&S) -> String,
    ) -> String {
        let mut r = String::new();

        r.push_str(&self.format_leaf_newick(i, f_leaf, f_tag));
        r.push_str(";\n");

        r
    }
}

impl<T: std::clone::Clone + std::fmt::Debug, S: std::fmt::Debug> PTree<T, S> {
    fn rec_disp(&self, i: usize, depth: usize) {
        println!(
            "{}{}#{:?} > CHILDREN: {:?}",
            " ".repeat(depth),
            i,
            self[i].tag,
            self[i].children
        );
        println!("{}    |- CONTENT: {:?}", " ".repeat(depth), self[i].content);

        for j in self[i].children.iter() {
            self.rec_disp(*j, depth + 2);
        }
    }
    pub fn disp(&self) {
        self.rec_disp(1, 0);
    }
}
