from pathlib import Path
import typing as T
import pandas as pd
from ..utils import get_logger
from .utils import parse_expression
from .expression import ExpressionError, check_probe_expressions, evaluate_expression
from collections import deque

log = get_logger(__name__)

class DAG():

    def __init__(self) -> None:
        self.nodes: T.List[Node] = []

    def from_config(self, config: dict):
        assert 'probes' in config, "probes key not found"
        assert isinstance(config['probes'], dict), "probes should be a dict"
        for probe_name, probe_config in config['probes'].items():
            assert 'template' in probe_config or 'expr' in probe_config, \
                "template or expr key not found"
            probe = TemplateProbe(self, probe_name, probe_config)
            self.nodes.append(probe)
        for node in self.nodes:
            if isinstance(node, ExprProbe):
                node.parse_expr()

    def get_downstream_nodes(self, node: "Node") -> T.List["Node"]:
        downstream = []
        for probe in self.nodes:
            if node in probe.deps:
                downstream.append(probe)
        return downstream

    def get_node_by_name(self, name: str) -> T.Optional["Node"]:
        for probe in self.nodes:
            if probe.name == name:
                return probe
        return None

    def get_all_0_deps(self) -> T.List["Node"]:
        return [probe for probe in self.nodes if len(probe.deps) == 0]

    def run(self, context: dict):
        """Run the DAG for a given context."""     
        for node in self.nodes:
            node.done = False
        queue = deque(self.get_all_0_deps())     
        built_count = 0
        while queue:
            node = queue.popleft()           
            all_deps_done = all(dep.done for dep in node.deps)           
            if not node.done and all_deps_done: 
                try:
                    node.build(context)
                    built_count += 1
                except Exception as e:
                    log.error(f"error building node {node.name}: {e}")
                    continue 
                for down_node in self.get_downstream_nodes(node):
                    if not down_node.done and all(dep.done for dep in down_node.deps):
                         if down_node not in queue:
                            queue.append(down_node)        


class Node:
    def __init__(self, dag: DAG, name: str, config: dict):
        self.dag = dag
        self.name = name
        self.config = config
        self.deps: T.List[Probe] = []
        self.done = False
        self.result = None 

    def build(self, context: dict):
        pass

class Probe(Node):
    pass

class ExprProbe(Probe):
    def __init__(self, dag: DAG, name: str, config: dict):
        super().__init__(dag, name, config)
        assert 'expr' in config, "expr key not found"
        self.expr = config['expr']
        self.external_deps: T.List[str] = []

    def parse_expr(self):
        try:
            deps = parse_expression(self.expr)
        except ExpressionError as e:
            raise ExpressionError(f"probe part {self.name!r}: {e}") from None
        for dep in deps:
            node = self.dag.get_node_by_name(dep)
            if node is None:
                self.external_deps.append(dep)
            else:
                self.deps.append(node)
    
    def build(self, context: dict):
        # Dependencies override context names; dotted part names such as
        # `probe_1.part1` are looked up as-is by the expression evaluator.
        namespace = context.copy()
        for dep in self.deps:
            if not dep.done:
                raise RuntimeError(f"dependency {dep.name} was not built before {self.name}")
            if dep.result is None:
                raise RuntimeError(f"error reading result for dependency {dep.name}: "
                                   f"dependency {dep.name} has no result in memory")
            namespace[dep.name] = str(dep.result)
        try:
            result = evaluate_expression(self.expr, namespace)
            self.result = str(result)  # Store only in memory
            self.done = True
            #log.info(f"successfully built {self.name}")
        except Exception as e:
            log.error(f"error evaluating {self.name}: {e}")
            raise

class TemplateProbe(Probe):
    def __init__(self, dag: DAG, name: str, config: dict):
        super().__init__(dag, name, config)
        assert 'template' in config, "template key not found"
        assert 'parts' in config, "parts key not found"
        self.parts: T.List[Probe] = []
        self.template: str = config['template']
        self.resolve_parts()

    def resolve_parts(self) -> None:
        self.parts.clear()
        for part_name, part_config in self.config['parts'].items():
            part: Node
            new_name = f"{self.name}.{part_name}"
            if 'template' in part_config:
                part = TemplateProbe(
                    self.dag, new_name, part_config)
            elif 'expr' in part_config:
                part = ExprProbe(
                    self.dag, new_name, part_config)
            else:
                raise ValueError(f"invalid part config: {part_config}")
            self.parts.append(part)
            self.dag.nodes.append(part)
        self.deps = self.parts

    def build(self, context: dict):
        if not all(part.done for part in self.parts):
             missing_deps = [p.name for p in self.parts if not p.done]
             raise RuntimeError(f"cannot build {self.name}, missing dependencies: {missing_deps}")
        part_results = {}
        for part in self.parts:
            part_name = part.name.split('.')[-1]
            if part.result is not None:
                part_results[part_name] = part.result
            else:
                raise RuntimeError(f"part {part.name} has no result in memory")
        try:
            result = self.template.format(**part_results)
            self.result = result  # Store only in memory
            self.done = True
        except Exception as e:
            log.error(f"error formatting template {self.name}: {e}")
            self.done = False
            raise

    def __getitem__(self, item: str):
        for part in self.parts:
            if part.name == item:
                return part

def construct_probes(config, contexts):
    """
    Construct probes for each target context.
    """
    # Reject invalid or unsafe expressions once, before any target is processed.
    errors = check_probe_expressions(config.get('probes'), check_names=False)
    if errors:
        raise ExpressionError("invalid probe expressions:\n" + "\n".join(errors))
    if not isinstance(contexts, list):
        context = contexts
        dag = DAG()
        dag.from_config(config)
        dag.run(context)
        return pd.DataFrame() 
    probe_results = []
    for context in contexts:
        try:
            dag = DAG()
            dag.from_config(config)
            dag.run(context)
            probe_data = {}
            for node in dag.nodes:
                if node.done and node.result is not None:
                    probe_data[node.name] = node.result
                else:
                    log.warning(f"node {node.name} not completed or result missing")
                    probe_data[node.name] = None
            probe_results.append(probe_data)
        except Exception as e:
            log.error(f"error processing target {context}: {e}")
            probe_data = {}
            if 'dag' in locals():
                probe_data = {node.name: None for node in dag.nodes}
            probe_results.append(probe_data)
    if probe_results:
        probe_df = pd.DataFrame(probe_results)
        probe_df = probe_df.fillna('')
    else:
        probe_df = pd.DataFrame()
    return probe_df 


