import copy


def string_composition(text, k):
    i = 0
    string_comp = []
    while i <= len(text) - k:
        string = text[i : i + k]
        string_comp.append(string)
        i = i + 1
    return string_comp


def string_spelled_by_genome_path(genome_path):
    k = len(genome_path[0])
    reconstruction = [genome_path[0][0 : k - 1]]
    for item in genome_path:
        reconstruction.append(item[k - 1])
    string_reconstruction = "".join(reconstruction)
    return string_reconstruction


def overlap_graph_problem(patterns):
    k = len(patterns[0])
    prefix_dict = {}
    adjacency_list = {}
    for item in patterns:
        key = item[0 : k - 1]
        val = prefix_dict.get(key)
        val = val + " " if val else ""
        val = val + item
        prefix_dict[key] = val
    print(prefix_dict)
    for key in patterns:
        suffix = key[1:k]
        if prefix_dict.get(suffix):
            adjacency_list[key] = prefix_dict.get(suffix)
    return adjacency_list


def de_bruijn_graph(patterns):
    de_bruijn_dict = {}
    k = len(patterns[0])
    for item in patterns:
        key = item[0 : k - 1]
        value = item[1:k]
        val = de_bruijn_dict.get(key)
        if val:
            val += " "
        else:
            val = ""
        val = val + value
        de_bruijn_dict[key] = val
    return de_bruijn_dict


def de_bruijn_graph_from_text(text, k):
    patterns = string_composition(text, k)
    print(patterns)
    de_bruijn_dict = de_bruijn_graph(patterns)
    return de_bruijn_dict


def eulerian_cycle(adjacency_dict: dict):
    current_cycle = []
    previous_cycle = []
    keys = [*adjacency_dict]
    next_node = keys[0]
    current_cycle.append(next_node)
    total_edges = 0
    for key in adjacency_dict:
        total_edges += len(adjacency_dict[key])

    while True:
        temp_dict = copy.deepcopy(adjacency_dict)
        keys = [*temp_dict]
        while True:
            outgoing_vertices = temp_dict[next_node]
            if len(outgoing_vertices) == 0 or len(current_cycle) == total_edges + 1:
                break
            from_previous_cycle = False
            for vertex in outgoing_vertices:
                if vertex in previous_cycle:
                    next_node = vertex
                    outgoing_vertices.remove(next_node)
                    from_previous_cycle = True
            if not from_previous_cycle:
                next_node = outgoing_vertices.pop(0)
            current_cycle.append(next_node)

        if len(current_cycle) == total_edges + 1:
            return current_cycle

        max_len = 0
        for key in current_cycle:
            value_array = temp_dict.get(key)
            if len(value_array) > max_len:
                max_len = len(temp_dict.get(key, []))
                next_node = key
        previous_cycle = copy.deepcopy(current_cycle)

    print("Current Cycle: ", current_cycle)
