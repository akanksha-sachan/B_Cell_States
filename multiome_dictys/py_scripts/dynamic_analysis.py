import matplotlib.pyplot as plt
import matplotlib.animation
from dictys.net import dynamic_network
from dictys.plot import layout, panel

def load_data(data_file):
    """
    Load the dynamic network data and define trajectory branches.
    
    Parameters:
        data_file (str): Path to the dynamic network file.
        branch_definitions (dict): Dictionary defining branches with start and end nodes.

    Returns:
        dynamic_network: Loaded dynamic network data.
        dict: Branches defining the network paths.
    """
    d0 = dynamic_network.from_file(data_file)
    return d0

def configure_branch_tfs(branch_name, tf_subnet, tf_annotations, target_annotations):
    """
    Configure the specific branch, transcription factors (TFs), and annotations for the dynamic subnetwork.

    Parameters:
        branch_name (str): Name of the branch being analyzed.
        tf_subnet (list of lists): TFs for each row's dynamic subnetwork graph.
        tf_annotations (list of lists): TFs for each row's other plots.
        target_annotations (list): Genes to annotate as targets in all rows.

    Returns:
        dict: Dictionary containing branch-specific configuration parameters.
    """
    config = {
        'branch_name': branch_name,
        'tfs_subnet': tf_subnet,
        'tfs_ann': tf_annotations,
        'target_ann': target_annotations
    }
    return config

def generate_and_save_animation(data, branch, config, vrange, nframe=100, fps=0.10, dpi=100, output_path='animation.mp4'):
    """
    Generate the dynamic network animation and save it to a file.

    Parameters:
        data: Loaded dynamic network data.
        branch (tuple): Start and end nodes defining the branch.
        config (dict): Configuration dictionary with TFs and target genes.
        vrange (dict): Value range for coloring in visualizations.
        nframe (int): Number of frames for the animation (default: 100).
        fps (float): Frames per second for the animation (default: 0.10).
        dpi (int): DPI for the animation (default: 100).
        output_path (str): Path to save the animation file.

    Returns:
        None
    """
    # Draw the dynamic network with notch layout
    layout1 = layout.notch(nframe=nframe, dpi=dpi)
    
    pts, fig, panels, animate_ka = layout1.draw(
        data, branch,
        bcde_tfs=config['tfs_ann'],
        e_targets=config['target_ann'],
        f_tfs=config['tfs_subnet'],
        a_ka={'scatterka': {'legend_loc': (0.6, 1)}},
        # Use the value range for coloring
        e_ka={'lim': vrange.get('Switching time', [-0.02, 0.02])},  # Example for 'Switching time'
    )
    
    ca = panel.animate_generic(pts, fig, panels)
    anim = ca.animate(**animate_ka)

    # Save the animation
    writer = matplotlib.animation.writers['ffmpeg_file'](fps=fps * nframe, codec='mpeg4')
    writer.frame_format = 'jpeg'
    anim.save(output_path, writer=writer, dpi='figure')

    print(f'Animation saved to {output_path}')


# Example usage of the functions
if __name__ == '__main__':
    # Load data
    data_file = '/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs/output/dynamic.h5'
    # Define branches with start and end nodes (2: activated b-cells, 1 is PB, 3 is GC)
    branches = {
        'Plasma-Blast': (2, 1),
        'Germinal-Center': (2, 3)
    }
    # Define value ranges for coloring
    vrange = {
        'Terminal logFC': [-4.5, 4.5],
        'Transient logFC': [-1.5, 1.5],
        'Switching time': [0.001, 0.01],
    }
    data = load_data(data_file)

    # Step 2: Configure the branch-specific transcription factors (TFs) and annotations
    branch_names = ['Plasma-Blast', 'Germinal-Center']
    tf_subnet = [['IRF4'], ['IRF8']]
    tf_annotations = [['IRF4', 'IRF8'], ['PRDM1', 'BCL6']]
    target_annotations = ['IRF4', 'IRF8', 'PRDM1', 'BCL6']
    config = configure_branch_tfs(branch_names[0], tf_subnet, tf_annotations, target_annotations)

    # Step 3: Generate and save the animation
    branch_pb = branches[branch_names[0]]
    branch_gc = branches[branch_names[1]]
    output_path_pb = f'/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs/output/animation-{branch_names[0]}.mp4'
    output_path_gc = f'/ocean/projects/cis240075p/asachan/datasets/B_Cell/multiome_1st_donor_UPMC_aggr/dictys_outs/output/animation-{branch_names[1]}.mp4'
    generate_and_save_animation(data, branch_pb, config, vrange, nframe=100, fps=0.10, dpi=100, output_path=output_path_pb)
    generate_and_save_animation(data, branch_gc, config, vrange, nframe=100, fps=0.10, dpi=100, output_path=output_path_gc)
