"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import pixel_shift
import pyomo.environ as pyomo
import geo_sub_priority_multi as geosub
import dill as pickle
import run_all_utilities as utils
import argparse

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


def main():
    print "\n"
    print "Copyright 2016 Sandia Corporation. Under the terms of Contract"
    print "DE-AC04-94AL85000, there is a non-exclusive license for use of"
    print "this work by or on behalf of the U.S. Government. Export of this"
    print "program may require a license from the United States Government."
    print "\n"

    parser = argparse.ArgumentParser(description='Specify output directory.')
    parser.add_argument('--o', dest='outdir', type=str, nargs='?', default='big_run_results',
                        help='an output directory, dude')
    parser.add_argument('-i', '--input_file', default="examples/mississippi/rivers160.bmp",
                        help='Specify the input bitmap image')
    parser.add_argument('-s', '--solver', default='cplex', help='Specify the solver for Pyomo to use')
    parser.add_argument('-bw', '--bandwidth', default=24, help='Specify the bandwidth constraint in terms of chips')
    parser.add_argument('-ns', '--num_subfoots', default=6, help='Specify the number of subfootprints for the problem')
    parser.add_argument('-ncr', '--num_chips_per_row', default=10,
                        help='Specify the number of chips in one row of an image')
    parser.add_argument('-ncc', '--num_chips_per_col', default=10,
                        help='Specify the number of chips in one column of an image')

    args = parser.parse_args()
    print(args.outdir)

    # Get the command line arguments
    outdir, bw, num_subfoots, filename, chip_shape, solver = ParseArguments(args)

    num_chips_per_row = chip_shape[0]
    num_chips_per_col = chip_shape[1]

    # Step 1: Read in the pixel list from a bitmap
    # read bitmap, its dimensions are pixel_shape, its list of non-zeros is pixel_list and those values are pixel_priorities
    pixel_list, pixel_priorities, pixel_shape = pixel_shift.BitmapToPixelList(filename)

    # Step 2 and 3: Get all possible chip layouts from shifted pixels positions.

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

    for index, chip_layout in enumerate(all_chip_layouts):
        m_title = args.outdir + "/MIP Results for Layout " + str(index)

        chip_priorities = pixel_shift.GenerateChipPriorities(chip_layout, pixel_configuration_list[index],
                                                             pixel_priorities, pixels_per_chip_row, pixels_per_chip_col)

        chips_to_cover = pixel_shift.ChipsToCover(chip_layout)

        print("Creating MIP model...")
        M = geosub.pyomo_create_model(None, None, chip_priorities, chips_to_cover, num_subfoots, bw, chip_shape[0],
                                      chip_shape[1])

        print("Creating MIP instance...")
        instance = M.create()

        print("Solving MIP instance...")
        results = opt.solve(instance)
        instance.load(results)
        m_title += "-Value: " + str(instance.obj_priority())
        geosub.pyomo_postprocess(None, instance, chip_priorities, chips_to_cover, num_subfoots, chip_layout,
                                 pixel_configuration_list[index], pixel_priorities, pixel_shape[0], pixel_shape[1],
                                 m_title)

        pickle_fn = args.outdir + "/mip_" + str(index) + ".pkl"
        output_pickle = open(pickle_fn, 'wb')
        pickle.dump(instance, output_pickle)
        pickle.dump(results, output_pickle)


if __name__ == "__main__":
    main()
