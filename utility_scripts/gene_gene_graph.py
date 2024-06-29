import pandas as pd

import pandas as pd
import networkx as nx

from tqdm import tqdm
AD_with_Plaques_df = pd.read_csv('preprocessed_all_stratification_grn/AD_with_Plaques_df.csv')
NCI_with_Plaques_df = pd.read_csv('preprocessed_all_stratification_grn/NCI_with_Plaques_df.csv')


for s in AD_with_Plaques_df.columns[2:]:
    B = nx.Graph()

    # Add nodes from the first column (e.g., genes)
    B.add_nodes_from(AD_with_Plaques_df['TF'])

    # Add nodes from the second column (e.g., TFs)
    B.add_nodes_from(AD_with_Plaques_df['Genes'])

    from tqdm import tqdm
    for index, row in tqdm(AD_with_Plaques_df.iterrows()):
        B.add_edge(row['TF'], row['Genes'], weight=row[s])
        
    edges = list(B.edges(data=True))

# Create a DataFrame from the list of edges
    df = pd.DataFrame(edges, columns=['Node1', 'Node2', 'Weight'])
    
    gene_nodes = list(set(df.Node2.values))
    df = pd.DataFrame()
    projected_graph = nx.MultiGraph()
    
   

# Create the projected graph with mean weights
    for u in tqdm(gene_nodes):
        for v in gene_nodes:
            if u != v:
                #print(u,v)
                u_neighbors = set(B[u])
                v_neighbors = set(B[v])
                common_neighbors = u_neighbors.intersection(v_neighbors)
                #print(common_neighbors)

                if len(common_neighbors) == 0:
                    continue
                #print(B[u][list(common_neighbors)[0]]['weight'])
                if common_neighbors:
                    common_neighbors = list(common_neighbors)
                    for i in common_neighbors:
                        mean_weight = sum(B[w][i]['weight'] for w in [u,v]) / 2
                        #print(mean_weight)
                        projected_graph.add_edge(u, v, weight=mean_weight)

    # Print the edges and their weights in the projected graph
    edges = list(projected_graph.edges(data=True))

# Create a DataFrame from the list of edges
    df = pd.DataFrame(edges, columns=['Node1', 'Node2', 'Weight'])
    
    df.to_csv('./Gene_Gene_Graph/AD_with_plaques/'+str(s)+ '.csv')
    break