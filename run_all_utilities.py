"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import perfect
import pixel_shift
import geo_sub_priority_multi as geosub
import heuristics as he

import dill as pickle

import os


def MakeResultsDictionary(chips_to_cover, chip_priorities, num_subfoots, chip_rows, chip_cols, chip_layout,
                           pixel_config, pixel_priorities, pixel_rows, pixel_cols):
    """
    Saves data into a dictionary for easy indexing and saving
    :param chips_to_cover:
    :param chip_priorities:
    :param num_subfoots:
    :param chip_rows:
    :param chip_cols:
    :param chip_layout:
    :param pixel_config:
    :param pixel_priorities:
    :return:
    """
    d = {"chips_to_cover": chips_to_cover, "chip_priorities": chip_priorities, "num_subfoots": num_subfoots, "chip_rows": chip_rows,
        "chip_cols": chip_cols, "chip_layout": chip_layout, "pixel_config": pixel_config, "pixel_priorities": pixel_priorities,
        "pixel_rows": pixel_rows, "pixel_cols": pixel_cols}
    return d


def MakeDirIfNotThere(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)


def GenerateSortedHeuristicLists(greedy_value_list, perfect_value_list):
    sorted_greedy_index_list = [i[0] for i in sorted(enumerate(greedy_value_list), key=lambda x: x[1])]
    sorted_perfect_index_list = [i[0] for i in sorted(enumerate(perfect_value_list), key=lambda x: x[1])]

    return sorted_greedy_index_list, sorted_perfect_index_list


def ComputeGreedyAndPerfectSolutions(all_chip_layouts, pixel_configuration_list, pixel_priorities, pixels_per_chip_row,
                                     pixels_per_chip_col, num_subfoots, bw, chip_shape):
    greedy_output = []
    greedy_value_list = []
    perfect_output = []
    perfect_value_list = []
    for index, chip_layout in enumerate(all_chip_layouts):

        print "chip layout %d of %d" % (index, len(all_chip_layouts) - 1)

        chip_priorities = pixel_shift.GenerateChipPriorities(chip_layout, pixel_configuration_list[index],
                                                             pixel_priorities,
                                                             pixels_per_chip_row, pixels_per_chip_col)

        chips_to_cover = pixel_shift.ChipsToCover(chip_layout)

        greedy_sol = he.greedy(num_subfoots, bw, chip_shape, chips_to_cover, chip_priorities)
        print greedy_sol
        perfect_sol = perfect.perfect(chip_priorities, bw)
        greedy_output.append(greedy_sol)
        perfect_output.append(perfect_sol)

        total_value = 0
        for sol in greedy_sol:
            total_value += sol[2]
        greedy_value_list.append(total_value)

        perfect_value_list.append(perfect_sol)

    return greedy_output, greedy_sol, greedy_value_list, perfect_output, perfect_sol, perfect_value_list


def PickleGreedyAndPerfect(greedy_output, greedy_value_list, perfect_output, perfect_value_list, sorted_index_list,
                           sorted_perfect_index_list, pickle_name):
    output_pickle = open(pickle_name, 'wb')
    pickle.dump(greedy_output, output_pickle)
    pickle.dump(greedy_value_list, output_pickle)
    pickle.dump(perfect_output, output_pickle)
    pickle.dump(perfect_value_list, output_pickle)
    pickle.dump(sorted_index_list, output_pickle)
    pickle.dump(sorted_perfect_index_list, output_pickle)


def SolveMIP(chip_priorities, chips_to_cover, num_subfoots, bw, chip_shape_row, chip_shape_col, solver):
    M = geosub.pyomo_create_model(None, None, chip_priorities, chips_to_cover, num_subfoots, bw, chip_shape_row,
                                  chip_shape_col)

    print("Creating MIP instance...")
    instance = M.create()

    print("Solving MIP instance...")
    results = solver.solve(instance)
    instance.load(results)
    return instance, results

def GenerateMipSolution(instance):
    """
    Generates the printable subfootprint solution for display.
    :param instance:
    :return: A list of subfootprints, each defined by its corners
    """
    subfeet = []
    for f in instance.F:
        ll_x = instance.Fkllx[f].value
        ll_y = instance.Fklly[f].value
        ur_x = instance.Fkurx[f].value
        ur_y = instance.Fkury[f].value
        subfootprint = [ll_x, ll_y, ur_x, ur_y]
        subfeet.append(subfootprint)
    return subfeet

def PostProcessMIP(instance, results, chip_priorities, chips_to_cover, num_subfoots, chip_layout, pixel_config,
                   pixel_priorities, pixel_shape_row, pixel_shape_col, plot_title, outdir, index):
    plot_title += "-Value: " + str(instance.obj_priority())

    # Plots MIP Results
    geosub.pyomo_postprocess(None, instance, chip_priorities, chips_to_cover, num_subfoots, chip_layout,
                             pixel_config,
                             pixel_priorities, pixel_shape_row, pixel_shape_col, plot_title)

    pickle_fn = outdir + "/mip_" + str(index) + ".pkl"
    output_pickle = open(pickle_fn, 'wb')
    pickle.dump(instance, output_pickle)
    pickle.dump(results, output_pickle)
