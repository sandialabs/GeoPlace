"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import dill as pickle
import pixel_shift
import pyomo.environ as pyomo
import post_process_multi as pp
import run_all_utilities as utils

"""
Our objectives are the following:
1: Read a bitmap image and get a corresponding pixel_list
2: Find all shifts of the pixel_list
3: Generate chip_layouts per each
4: Loop over the chip layouts, applying either the MIP or a Heuristic or some other solver
"""


def ParseArguments(args):
    outdir = args.outdir

    # Make the directory to store output data if it doesn't exist yet
    utils.MakeDirIfNotThere(args.outdir)
    bw = args.bandwidth
    num_subfoots = args.num_subfoots
    filename = args.input_file
    solver = args.solver

    # define number of chips in the footprint, in each direction
    # these are assumed to divide evenly into the image shape
    # e.g. 15 x 20 pixels = 3 x 4 chips of 5 x 5 pixels each

    num_chips_per_row = args.num_chips_per_row
    num_chips_per_col = args.num_chips_per_col
    chip_shape = (num_chips_per_row, num_chips_per_col)

    print("Input file is specified to be: ", filename)
    print("The bandwidth is specified to be", bw)
    print("The number of subfootprints is specified to be", num_subfoots)

    # each subfoot can have up to bw/num_subfoots chips
    print "each subfoot can have", bw / num_subfoots, "chips"
    print("The output directory is specified to be", outdir)
    return outdir, bw, num_subfoots, filename, chip_shape, solver


def GeneratePlotTitle(count, index, greedy_value_list, perfect_value_list, outdir):
    if count < 11:
        h_title = outdir + "/Greedy " + str(count) + " Results Layout " + str(index) + "-Value: " + str(
            greedy_value_list[index])
        m_title = outdir + "/MIP Results for Greedy Layout " + str(index)
    else:
        h_title = outdir + "/Perfect " + str(count - 10) + " Results Layout " + str(index) + "-Value: " + str(
            perfect_value_list[index])
        m_title = outdir + "/MIP Results for Perfect Layout " + str(index)
    return h_title, m_title


def main():
    print "\n"
    print "Copyright 2016 Sandia Corporation. Under the terms of Contract"
    print "DE-AC04-94AL85000, there is a non-exclusive license for use of"
    print "this work by or on behalf of the U.S. Government. Export of this"
    print "program may require a license from the United States Government."
    print "\n"

    import argparse

    parser = argparse.ArgumentParser(description='Specify output directory.')
    parser.add_argument('--o', dest='outdir', type=str, nargs='?', default='big_run_results',
                        help='an output directory, dude')
    parser.add_argument('-i', '--input_file', default="examples/mississippi/rivers160.bmp",
                        help='Specify the input bitmap image')
    parser.add_argument('-bw', '--bandwidth', default=24, help='Specify the bandwidth constraint in terms of chips')
    parser.add_argument('-ns', '--num_subfoots', default=6, help='Specify the number of subfootprints for the problem')
    parser.add_argument('-ncr', '--num_chips_per_row', default=10,
                        help='Specify the number of chips in one row of an image')
    parser.add_argument('-ncc', '--num_chips_per_col', default=10,
                        help='Specify the number of chips in one column of an image')
    parser.add_argument('-s', '--solver', default='cplex', help='Specify the solver for Pyomo to use')
    args = parser.parse_args()
    outdir, bw, num_subfoots, filename, chip_shape, solver = ParseArguments(args)

    # Step 1: Read in the pixel list from a bitmap
    # read bitmap, its dimensions are pixel_shape, its list of non-zeros is pixel_list and those values are pixel_priorities
    pixel_list, pixel_priorities, pixel_shape = pixel_shift.BitmapToPixelList(filename)

    # Step 2 and 3: Get all possible chip layouts from shifted pixels positions.

    num_chips_per_row = chip_shape[0]
    num_chips_per_col = chip_shape[1]
    pixels_per_chip_row = pixel_shape[0] / chip_shape[0]
    pixels_per_chip_col = pixel_shape[1] / chip_shape[1]

    print "domain is", pixel_shape[0], "x", pixel_shape[1], "pixels"
    print "footprint is", num_chips_per_row * pixels_per_chip_row, "x", num_chips_per_col * pixels_per_chip_col, "pixels"
    print "footprint is", num_chips_per_row, "x", num_chips_per_col, "chips"
    print "chip is", pixels_per_chip_row, "x", pixels_per_chip_col, "pixels"

    all_chip_layouts, pixel_configuration_list = pixel_shift.FindAllChipLayouts(pixel_list, pixels_per_chip_row,
                                                                                pixels_per_chip_col, chip_shape[0],
                                                                                chip_shape[1])

    opt = pyomo.SolverFactory(solver)

    ## Step 4: Loop over chip layouts

    print "Running greedy heuristic for pixel shift"

    greedy_output, greedy_sol, greedy_value_list, perfect_output, perfect_sol, perfect_value_list = utils.ComputeGreedyAndPerfectSolutions(
        all_chip_layouts, pixel_configuration_list, pixel_priorities, pixels_per_chip_row,
        pixels_per_chip_col, num_subfoots, bw, chip_shape)

    sorted_greedy_index_list, sorted_perfect_index_list = utils.GenerateSortedHeuristicLists(greedy_value_list,
                                                                                             perfect_value_list)

    pickle_name = args.outdir + "/greedies.pkl"
    utils.PickleGreedyAndPerfect(greedy_output, greedy_value_list, perfect_output, perfect_value_list,
                                 sorted_greedy_index_list, sorted_perfect_index_list, pickle_name)

    # Get top 10 best chip layouts from greedy values
    heuristic_chip_layouts = sorted_greedy_index_list[-1:]

    # Add on the top ten perfect values
    for val in sorted_perfect_index_list[-1:]:
        heuristic_chip_layouts.append(val)

    count = 0
    for index2 in heuristic_chip_layouts:
        count += 1
        chip_layout = all_chip_layouts[index2]
        pixel_config = pixel_configuration_list[index2]

        h_title, m_title = GeneratePlotTitle(count, index2, greedy_value_list, perfect_value_list, args.outdir)
        chip_priorities = pixel_shift.GenerateChipPriorities(chip_layout, pixel_config, pixel_priorities,
                                                             pixels_per_chip_row, pixels_per_chip_col)

        chips_to_cover = pixel_shift.ChipsToCover(chip_layout)

        greedy_sol = greedy_output[index2]

        # Plot the Greedy solution results
        pp.PlotChipResults(chips_to_cover, chip_priorities, num_subfoots, chip_shape[0], chip_shape[1], chip_layout,
                           pixel_config, pixel_priorities, pixel_shape[0], pixel_shape[1], h_title, pp.GreedySolutionConvert(greedy_sol))

        print greedy_output[index2]

        instance, results = utils.SolveMIP(chip_priorities, chips_to_cover, num_subfoots, bw, chip_shape[0],
                                           chip_shape[1], opt)

        mip_solution = utils.GenerateMipSolution(instance)

	m_title  += "-Value: " + str(instance.obj_priority())

        pp.PlotChipResults(chips_to_cover, chip_priorities, num_subfoots, chip_shape[0], chip_shape[1], chip_layout,
                   pixel_config, pixel_priorities, pixel_shape[0], pixel_shape[1], m_title,
                   mip_solution)

if __name__ == "__main__":
    main()
