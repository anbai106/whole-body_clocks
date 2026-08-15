import numpy as np
import pandas as pd
import pyvista as pv
import os
from matplotlib import pyplot as plt

def plot_in_pyvsita_with_bg(mesh_stat, mesh_bg_left, mesh_bg_right, min_stat, max_stat, output_name):
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
        fmt="%.10f",
        font_family="arial",
    )

    cpos_1 = [(0.6, -17, 514.9), (0.61, -17.1, 11.86), (0, 1, 0)] ### for axial view
    cpos_2 = [(-11.21, 451.62, -101.12), (1.78, -12.13, 0.86), (-0.003, 0.21, 0.97)]
    cpos_3 = [(1.21, 1.62, 101.12), (1.78, 1.13, 0.86), (0.003, 0.1, 0.97)]

    pv.set_plot_theme("document")
    plotter = pv.Plotter()
    plotter.add_mesh(mesh_stat, scalar_bar_args=sargs, clim=[min_stat, max_stat], cmap="coolwarm")
    plotter.add_mesh(mesh_bg_left, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter.add_mesh(mesh_bg_right, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter.show(cpos=cpos_1, screenshot=os.path.join(output_name + '_view1_top.png'))
    plotter.close(render=False)

    plotter2 = pv.Plotter()
    plotter2.add_mesh(mesh_stat, scalar_bar_args=sargs, clim=[min_stat, max_stat], cmap="coolwarm")
    plotter2.add_mesh(mesh_bg_left, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter2.add_mesh(mesh_bg_right, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter2.show(cpos=cpos_2, screenshot=os.path.join(output_name + '_view2_top.png'))
    plotter2.close(render=False)

    plotter3 = pv.Plotter()
    plotter3.add_mesh(mesh_stat, scalar_bar_args=sargs, clim=[min_stat, max_stat], cmap="coolwarm")
    plotter3.add_mesh(mesh_bg_left, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter3.add_mesh(mesh_bg_right, opacity=0.1, cmap="coolwarm")  # add a mesh to the scene
    plotter3.show(cpos=cpos_3, screenshot=os.path.join(output_name + '_view3_top.png'))
    plotter3.close(render=False)

def plot_in_pyvsita_only_bg(mesh_bg_left, mesh_bg_right, output_name):
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
    cpos_1 = [(0.6, -17, 514.9), (0.61, -17.1, 11.86), (0, 1, 0)] ### for axial view
    cpos_2 = [(-11.21, 451.62, -101.12), (1.78, -12.13, 0.86), (-0.003, 0.21, 0.97)]

    plotter = pv.Plotter()
    plotter.add_mesh(mesh_bg_left, opacity=0.1)  # add a mesh to the scene
    plotter.add_mesh(mesh_bg_right, opacity=0.1)  # add a mesh to the scene
    plotter.show(cpos=cpos_1, screenshot=os.path.join(output_name + '_view1_top.png'))
    plotter.close(render=False)

    plotter2 = pv.Plotter()
    plotter2.add_mesh(mesh_bg_left, opacity=0.1)  # add a mesh to the scene
    plotter2.add_mesh(mesh_bg_right, opacity=0.1)  # add a mesh to the scene
    plotter2.show(cpos=cpos_2, screenshot=os.path.join(output_name + '_view2_top.png'))
    plotter2.close(render=False)


def visual_stat_image(input_tsv, dictionary_tsv, JHU_label_tract_dir, bg_left, bg_right, output_dir):
    """
    This is to write the stat of each components into the NIFTI image for visualization
    Args:
        input_tsv:
        mist_atlas:

    Returns:

    """
    ## background meshes
    mesh_bg_left = pv.read(bg_left)
    mesh_bg_right = pv.read(bg_right)
    df_dict = pd.read_csv(dictionary_tsv, sep='\t')

    ### manipulate the intput file and merge them
    df_stat = pd.read_csv(input_tsv, sep='\t')
    df_stat = df_stat.loc[df_stat['IDP'].isin([
    "mean_fa_in_middle_cerebellar_peduncle_on_fa_skeleton_f25056_2_0",
    "mean_fa_in_pontine_crossing_tract_on_fa_skeleton_f25057_2_0",
    "mean_fa_in_genu_of_corpus_callosum_on_fa_skeleton_f25058_2_0",
    "mean_fa_in_body_of_corpus_callosum_on_fa_skeleton_f25059_2_0",
    "mean_fa_in_splenium_of_corpus_callosum_on_fa_skeleton_f25060_2_0",
    "mean_fa_in_fornix_on_fa_skeleton_f25061_2_0",
    "mean_fa_in_corticospinal_tract_on_fa_skeleton_right_f25062_2_0",
    "mean_fa_in_corticospinal_tract_on_fa_skeleton_left_f25063_2_0",
    "mean_fa_in_medial_lemniscus_on_fa_skeleton_right_f25064_2_0",
    "mean_fa_in_medial_lemniscus_on_fa_skeleton_left_f25065_2_0",
    "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25066_2_0",
    "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25067_2_0",
    "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25068_2_0",
    "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25069_2_0",
    "mean_fa_in_cerebral_peduncle_on_fa_skeleton_right_f25070_2_0",
    "mean_fa_in_cerebral_peduncle_on_fa_skeleton_left_f25071_2_0",
    "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25072_2_0",
    "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25073_2_0",
    "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25074_2_0",
    "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25075_2_0",
    "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25076_2_0",
    "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25077_2_0",
    "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_right_f25078_2_0",
    "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_left_f25079_2_0",
    "mean_fa_in_superior_corona_radiata_on_fa_skeleton_right_f25080_2_0",
    "mean_fa_in_superior_corona_radiata_on_fa_skeleton_left_f25081_2_0",
    "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_right_f25082_2_0",
    "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_left_f25083_2_0",
    "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25084_2_0",
    "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25085_2_0",
    "mean_fa_in_sagittal_stratum_on_fa_skeleton_right_f25086_2_0",
    "mean_fa_in_sagittal_stratum_on_fa_skeleton_left_f25087_2_0",
    "mean_fa_in_external_capsule_on_fa_skeleton_right_f25088_2_0",
    "mean_fa_in_external_capsule_on_fa_skeleton_left_f25089_2_0",
    "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25090_2_0",
    "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25091_2_0",
    "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_right_f25092_2_0",
    "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_left_f25093_2_0",
    "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25094_2_0",
    "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25095_2_0",
    "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25096_2_0",
    "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25097_2_0",
    "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25098_2_0",
    "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25099_2_0",
    "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_right_f25100_2_0",
    "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_left_f25101_2_0",
    "mean_fa_in_tapetum_on_fa_skeleton_right_f25102_2_0",
    "mean_fa_in_tapetum_on_fa_skeleton_left_f25103_2_0"])]
    # Define tract list in the correct order
    tracts = [
        "mean_fa_in_middle_cerebellar_peduncle_on_fa_skeleton_f25056_2_0",
        "mean_fa_in_pontine_crossing_tract_on_fa_skeleton_f25057_2_0",
        "mean_fa_in_genu_of_corpus_callosum_on_fa_skeleton_f25058_2_0",
        "mean_fa_in_body_of_corpus_callosum_on_fa_skeleton_f25059_2_0",
        "mean_fa_in_splenium_of_corpus_callosum_on_fa_skeleton_f25060_2_0",
        "mean_fa_in_fornix_on_fa_skeleton_f25061_2_0",
        "mean_fa_in_corticospinal_tract_on_fa_skeleton_right_f25062_2_0",
        "mean_fa_in_corticospinal_tract_on_fa_skeleton_left_f25063_2_0",
        "mean_fa_in_medial_lemniscus_on_fa_skeleton_right_f25064_2_0",
        "mean_fa_in_medial_lemniscus_on_fa_skeleton_left_f25065_2_0",
        "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25066_2_0",
        "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25067_2_0",
        "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25068_2_0",
        "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25069_2_0",
        "mean_fa_in_cerebral_peduncle_on_fa_skeleton_right_f25070_2_0",
        "mean_fa_in_cerebral_peduncle_on_fa_skeleton_left_f25071_2_0",
        "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25072_2_0",
        "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25073_2_0",
        "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25074_2_0",
        "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25075_2_0",
        "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25076_2_0",
        "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25077_2_0",
        "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_right_f25078_2_0",
        "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_left_f25079_2_0",
        "mean_fa_in_superior_corona_radiata_on_fa_skeleton_right_f25080_2_0",
        "mean_fa_in_superior_corona_radiata_on_fa_skeleton_left_f25081_2_0",
        "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_right_f25082_2_0",
        "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_left_f25083_2_0",
        "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25084_2_0",
        "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25085_2_0",
        "mean_fa_in_sagittal_stratum_on_fa_skeleton_right_f25086_2_0",
        "mean_fa_in_sagittal_stratum_on_fa_skeleton_left_f25087_2_0",
        "mean_fa_in_external_capsule_on_fa_skeleton_right_f25088_2_0",
        "mean_fa_in_external_capsule_on_fa_skeleton_left_f25089_2_0",
        "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25090_2_0",
        "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25091_2_0",
        "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_right_f25092_2_0",
        "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_left_f25093_2_0",
        "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25094_2_0",
        "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25095_2_0",
        "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25096_2_0",
        "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25097_2_0",
        "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25098_2_0",
        "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25099_2_0",
        "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_right_f25100_2_0",
        "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_left_f25101_2_0",
        "mean_fa_in_tapetum_on_fa_skeleton_right_f25102_2_0",
        "mean_fa_in_tapetum_on_fa_skeleton_left_f25103_2_0"
    ]
    # Make mapping DataFrame: tract name -> tract number
    tract_map = pd.DataFrame({
        "IDP": tracts,
        "Tract_number": list(range(1, len(tracts) + 1))
    })
    # Merge with your stats
    df_stat = df_stat.merge(tract_map, on="IDP", how="left")
    df_stat_OR = df_dict.merge(df_stat, how='inner', left_on='Tract_number', right_on='Tract_number')
    p_thres = 0.05/720/9

    ## Read stats
    group_list = ['OR']  ## for the abstract
    for metric in group_list:
        print("Visualize stats for %s" % metric)
        df = df_stat_OR.loc[df_stat_OR['P_value'] < p_thres]
        df = df.loc[df['OR'] > 1]
        ### select only the brain MRIBAG
        df = df[df['BAG'].isin(['Brain_MRIBAG'])]

        min_stat = df['OR'].min()
        max_stat = df['OR'].max()

        df_new = df.copy(deep=True)
        tract_list = list(list(df_new['Tract_abb']))

        j = 0
        for i in tract_list:
            df_row = df_new[df_new['Tract_abb'].isin([i])]
            tract_num = df_row.Tract_number.values[0]
            tract_abb = df_row.Tract_abb.values[0]
            stat_OR = df_row['OR'].values[0]
            print("Tract: %s is significant, OR-value: %f" % (tract_abb, stat_OR))
            vtk_file = os.path.join(JHU_label_tract_dir, 'JHU_' + tract_abb + '_' + str(tract_num) + '.vtk')
            if not os.path.exists(vtk_file):
                raise Exception("Sth wrong is here...")
            if j == 0:
                j += 1
                mesh_final = pv.read(vtk_file)
                mesh_final.cell_data['OR'] = np.ones(mesh_final.n_cells) * stat_OR
            else:
                j += 1
                mesh = pv.read(vtk_file)
                mesh.cell_data['OR'] = np.ones(mesh.n_cells) * stat_OR
                ### concatenate to the finla mesh_final
                mesh_final += mesh
        ### plot the figure
        if j != 0:
            mesh_final.save(os.path.join(JHU_label_tract_dir, 'JHU_Label_' + metric + '_pyvista.vtk'))
        if j != 0:
            output_name = os.path.join(output_dir, 'JHU_Label_' + metric)
            plot_in_pyvsita_with_bg(mesh_final, mesh_bg_left, mesh_bg_right, min_stat, max_stat, output_name)
            del mesh_final
        else:
            output_name = os.path.join(output_dir, 'JHU_Label_' + metric)
            plot_in_pyvsita_only_bg(mesh_bg_left, mesh_bg_right, output_name)

input_tsv = "/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Brain_WM.tsv"
dictionary_tsv = "/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/dictionary_JHU_labels.tsv"
JHU_label_tract_dir = "/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/atlas-JHULabel"
bg_left = '/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/FsAverage_hemi-left_pial.vtk'
bg_right = '/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/FsAverage_hemi-right_pial.vtk'
output_dir = "/Users/hao/Dropbox/2025_RSS/Fig/Orig"
visual_stat_image(input_tsv, dictionary_tsv, JHU_label_tract_dir, bg_left, bg_right, output_dir)