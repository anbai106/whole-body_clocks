import numpy as np
import pandas as pd
import pyvista as pv
import os
import re
def plot_in_pyvsita_with_bg(mesh_stat, mesh_bg, min_stat, max_stat, output_name):
    """
    Plots for top N components
    Args:
        mesh_stat:
        mesh_bg_left:
        mesh_bg_right:

    Returns:

    """
    sargs = dict(
        title_font_size=20,
        label_font_size=16,
        shadow=True,
        n_labels=3,
        italic=True,
        fmt="%.1f",
        font_family="arial",
    )
    pv.set_plot_theme("document")
    plotter = pv.Plotter()
    plotter.add_mesh(mesh_stat, scalar_bar_args=sargs, clim=[min_stat, max_stat], cmap="coolwarm")
    plotter.add_mesh(mesh_bg, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene

    ### for MUSE
    cpos_1 = [(371, -101, 99.9), (1.78, -12.1, 0.86), (0.20, 0.97, 0.09)] ### for temporal lobe
    cpos_2 = [(-11.21, 451.62, -101.12), (1.78, -12.13, 0.86), (-0.003, 0.21, 0.97)] ### for subcortical structures

    ### for AAL
    # cpos_1 = [(-400.88, -126.30, -112.34), (1.48, -9.02, -0.13), (-0.22, 0.95, -0.13)] ### for temporal lobe
    # cpos_2 = [(1.65, 420.15, 63.57), (1.48, -9.02, -0.13), (0.008, 0.146, -0.99)] ### for subcortical structures

    # def my_cpos_callback():
    #     plotter.add_text(str(plotter.camera_position), name="cpos")
    #     return
    # plotter.add_key_event("p", my_cpos_callback)
    plotter.show(cpos=cpos_1, screenshot=os.path.join(output_name + '_view1_all.png'))
    plotter.close(render=False)

    plotter2 = pv.Plotter()
    plotter2.add_mesh(mesh_stat, scalar_bar_args=sargs, clim=[min_stat, max_stat], cmap="coolwarm")
    plotter2.add_mesh(mesh_bg, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter2.show(cpos=cpos_2, screenshot=os.path.join(output_name + '_view2_all.png'))
    plotter2.close(render=False)

def visual_stat_image_aal2(input_tsv, muse_vtk_dir, C_string, bg, output_dir):
    """
    This is to write the stat of each components into the NIFTI image for visualization
    Args:
        input_tsv:
        MUSE_atlas:

    Returns:

    """
    ## background meshes
    mesh_bg = pv.read(bg)

    print("MUSE atlas for %s" % C_string)
    ## Read stats
    task_list = ['OR']
    df_stat = pd.read_csv(input_tsv, sep='\t')
    for task in task_list:
        print("%s for task %s" % (C_string, task))
        df_new = df_stat.copy(deep=True)
        df_new_sig = df_new[df_new['P_value'] < 0.05/720/9]
        df_new_sig = df_new_sig.loc[df_new_sig['OR'] > 1]
        ### select only the brain MRIBAG
        df_new_sig = df_new_sig[df_new_sig['BAG'].isin(['Brain_MRIBAG'])]
        min_stat = df_new_sig[task].min()
        max_stat = df_new_sig[task].max()
        ROI_list = list(df_new_sig.IDP)

        ### Read each vtk that was created by the Freesurfer, which has only geometry and topology info, Points, but not attributes, scalar for stats values
        if os.path.exists(os.path.join(muse_vtk_dir, C_string + '_' + task + '_pyvista.vtk')):
            mesh_final = pv.read(os.path.join(muse_vtk_dir, C_string + '_' + task + '_pyvista.vtk'))
        else:
            for i in range(df_new_sig.shape[0]):
                ROI_name = ROI_list[i]
                vtk_file = os.path.join(muse_vtk_dir, re.sub(r'^MUSE_Volume_(\d+)', r'MUSE_GM_\1', ROI_name) + '.vtk')
                if ROI_name in list(ROI_list):
                    row = df_new_sig.loc[df_new_sig['IDP'].isin([ROI_name])]
                    stat = row[task].values[0]
                    if i == 0:
                        mesh_final = pv.read(vtk_file)
                        mesh_final.cell_data['OR'] = np.ones(mesh_final.n_cells) * stat
                    else:
                        mesh = pv.read(vtk_file)
                        mesh.cell_data['OR'] = np.ones(mesh.n_cells) * stat
                        ### concatenate to the finla mesh_final
                        mesh_final += mesh
                else:
                    continue

            mesh_final.save(os.path.join(muse_vtk_dir, C_string + '_' + task + '_pyvista.vtk'))
        output_name = os.path.join(output_dir, task + '_' + C_string)
        plot_in_pyvsita_with_bg(mesh_final, mesh_bg, min_stat, max_stat, output_name)

input_tsv = "/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Brain_GM.tsv"
muse_vtk_dir = "/Users/hao/cubic-home/Dataset/Atlas/MUSE/VTK_GM"
C_string = 'MUSE'
bg = '/Users/hao/cubic-home/Dataset/Atlas/MUSE/VTK_GM/MUSE_GM_all.vtk'
output_dir = "/Users/hao/Dropbox/2025_RSS/Fig/Orig"
visual_stat_image_aal2(input_tsv, muse_vtk_dir, C_string, bg, output_dir)