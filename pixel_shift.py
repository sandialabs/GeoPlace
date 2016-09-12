"""
Copyright 2016 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000, there is a non-exclusive license for use of
this work by or on behalf of the U.S. Government. Export of this
program may require a license from the United States Government.
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm
import numpy as np
import PIL.Image as Image

"""
The purpose of this file is to take in a list of pixels of interest and compute all possible layouts of corresponding
chips.

We use the following convention: A pixel is represented as a (row,col). So, if pixel = (4,19), it is in row 4, column 19.
The (row,col) convention follows matrix conventions, so row 0 is the top row, and row N would be the bottom row.
Column 0 is the leftmost column, and column M would be M columns to the right. This allows us to use bitmaps and
numpy without playing any games with reindexing.

Chips will follow the same convention. If we say we have a (150,1000) array of chips, that means there are 150 rows with
1000 columns.
"""

def PixelInChip(pixel_row, pixel_col, num_rows_of_chips, num_cols_of_chips):
    """

    :param pixel_row: row-coordinate of activated pixel
    :param pixel_col: col-coordiante of activated pixel
    :param pixels_in_one_row_of_chip: number of pixels in a single row of a chip
    :param pixels_in_one_col_of_chip: number of pixels in a single column of a chip
    :return: A pair, chip_x, chip_y showing the x and y coordinates of the chip the pixel is in
    """
    # Given the pixel positions, which chip is the pixel in?

    # If there are 5 chips total, then floor(pixel_x / 5) shows which chip position the pixel is in
    chip_row = np.floor(pixel_row / num_rows_of_chips)

    chip_col = np.floor(pixel_col / num_cols_of_chips)

    return chip_row, chip_col


def FindActivatedChips(pixel_list, pixels_in_chip_row_shape, pixels_in_chip_col_shape, num_rows_of_chips, num_cols_of_chips):
    """

    :param pixel_list: A list of tuples of activated pixels
    :return: A list of which chips contain pixels
    """

    chip_list = np.zeros((num_rows_of_chips, num_cols_of_chips))
    for pixel in pixel_list:
        # Which chip is this pixel in?
        # A chip is shape (pixel_rows, pixel_cols)
        chip_row, chip_col = PixelInChip(pixel[0], pixel[1], pixels_in_chip_row_shape, pixels_in_chip_col_shape)
        chip_list[chip_row][chip_col] = 1  # Indicate this is an activated chip

    return chip_list


def ShiftPixels(pixel_list, shift_row, shift_col):
    """
    The purpose of this function is to shift all the pixels in a direction. If a pixel goes off of the FPA,
    it should fail. This case should be prevented by another method keeping track of the pixels.
    :param pixel_list: List of incoming pixels. Pixel (4,9) indicates row 4, column 9
    :param shift_row: How much to shift the row by. This is a vertical shift.
    :param shift_col: How much to shift the column by. This is a horizontal shift.
    :return: A shifted list of pixels
    """
    shifted_pixel_list = []
    for pixel in pixel_list:
        pixel_x = pixel[0] + shift_row
        pixel_y = pixel[1] + shift_col
        shifted_pixel_list.append((pixel_x, pixel_y))

    return shifted_pixel_list


def FindMaxShifts(pixel_list, rows_of_pixels, cols_of_pixels):
    """
    The purpose of this function is, given an initial configuration of activated pixels, what is the most we can shift
    to the left, to the right, up, and down? At some point, if we shift too many pixels, e.g., left, then at least one
    pixel will no longer be on the FPA, and hence we have lost information
    :param pixel_list: A list of pixels of interest in some initial configuration
    :param total_pixels_x: Total number of pixels in the x direction
    :param total_pixels_y: Total number of pixels in the y direction
    :return: The total number of left shifts, right shifts, down shifts, and up shifts we may do to the initial pixel
    list without losing any pixels off the sides of the FPA
    """

    # Cast as a numpy array
    pixel_array = np.array(pixel_list)
    min_pixel_row = np.min(pixel_array[:,0])
    max_pixel_row = np.max(pixel_array[:,0])

    min_pixel_col = np.min(pixel_array[:,1])
    max_pixel_col = np.max(pixel_array[:,1])
    print"Non-zero pixel row range is", min_pixel_row, "to", max_pixel_row, "out of 0 to",rows_of_pixels-1
    print"Non-zero pixel col range is", min_pixel_col, "to", max_pixel_col, "out of 0 to",cols_of_pixels-1
    left_shift = min_pixel_col
    print"The array of chips span", rows_of_pixels, " x ", cols_of_pixels, " pixels."
    right_shift = cols_of_pixels - max_pixel_col - 1
    down_shift = rows_of_pixels - max_pixel_row -1
    up_shift = min_pixel_row

    return left_shift, right_shift, down_shift, up_shift


def FindAllPixelConfigurations(pixel_list, left_shift_max, right_shift_max, down_shift_max, up_shift_max):
    """
    Given an initial list of pixels, and the maximal possible shifts (from initial configuration), compute all different
    shifted pixel configurations.
    :param pixel_list:
    :param left_shift_max: The furthest left we can shift before a pixel goes off the FPA (from initial configuration)
    :param right_shift_max: The furthest right we can shift before a pixel goes off the FPA (from initial configuration)
    :param down_shift_max: The furthest down we can shift before a pixel goes off the FPA (from initial configuration)
    :param up_shift_max: The furthest up we can shift before a pixel goes off the FPA (from initial configuration)
    :return:
    """

    print("Left", left_shift_max, "Right", right_shift_max, "Down", down_shift_max , "Up", up_shift_max)
    pixel_configuration_list = []
    for up_shift in range(0, up_shift_max+1):
        print("Going up by", up_shift)

        # Pixels should be shifted vertically first, horizontally second
        # TODO: NOTE THIS IS CONFUSING. A "positive" vertical shift is DOWN. For example, row 5 is beneath row 4.
        for left_shift in range(0, left_shift_max+1):
            new_pixel_list = ShiftPixels(pixel_list, -1*up_shift, -1 * left_shift)
            pixel_configuration_list.append(new_pixel_list)

        for right_shift in range(1, right_shift_max+1):
            new_pixel_list = ShiftPixels(pixel_list, -1*up_shift, right_shift)
            pixel_configuration_list.append(new_pixel_list)

    for down_shift in range(1, down_shift_max+1):
        print("Going down by", down_shift)

        for left_shift in range(0, left_shift_max+1):
            new_pixel_list = ShiftPixels(pixel_list,  down_shift, -1 * left_shift)
            pixel_configuration_list.append(new_pixel_list)

        for right_shift in range(1, right_shift_max+1):
            new_pixel_list = ShiftPixels(pixel_list,  down_shift, right_shift)
            pixel_configuration_list.append(new_pixel_list)

    return pixel_configuration_list


def FindAllChipLayouts(pixel_list, pixels_per_chip_row, pixels_per_chip_col, num_chips_row, num_chips_col):
    """
    From an initial list of pixels of interest, compute every feasible layout of chips by shifting the pixels in all
    possible directions.
    :param pixel_list: A list of activated pixels from the full image
    :param pixels_per_chip_row: The shape, in pixels, of a single chip in rows
    :param pixels_per_chip_col: The shape, in pixels, of a single chip in cols
    :param num_chips_row: The shape of the chip mask of the full image in rows
    :param num_chips_col: The shape of the chip mask of the full image in cols
    :return: A layout of chips
    """
    total_pixels_x = pixels_per_chip_row * num_chips_row
    total_pixels_y = pixels_per_chip_col * num_chips_col

    left_shift, right_shift, down_shift, up_shift = FindMaxShifts(pixel_list, total_pixels_x, total_pixels_y)

    # samitch: limit shifts to one period
    right_shift= min( right_shift, pixels_per_chip_col-1 )
    left_shift = min(  left_shift, pixels_per_chip_col-1 - right_shift)
    up_shift   = min(    up_shift, pixels_per_chip_row-1 )
    down_shift = min(  down_shift, pixels_per_chip_row-1 - up_shift)

    pixel_configuration_list = FindAllPixelConfigurations(pixel_list, left_shift, right_shift, down_shift, up_shift)

    chip_layout_list = []

    for pixel_list in pixel_configuration_list:
        activated_chips = FindActivatedChips(pixel_list, pixels_per_chip_row, pixels_per_chip_col, num_chips_row, num_chips_col)
        chip_layout_list.append(activated_chips)

    return chip_layout_list, pixel_configuration_list

def GenerateChipPriorities(chip_list, pixel_list, pixel_priorities, pixels_per_chip_row_shape, pixels_per_chip_col_shape):
    """

    :param chip_list: A numpy array of shape (rows,cols) indicating whether a chip is activated or not
    :param pixel_list: List of pixels of interest
    :param pixel_priorities: List of pixel priorities corresponding to pixels in pixel_list
    :return: A list of priorities for each chip
    """


    chip_priorities = np.zeros(np.shape(chip_list))
    for index, pixel in enumerate(pixel_list):
        # Which chip is it in?
        chip_row, chip_col = PixelInChip(pixel[0], pixel[1], pixels_per_chip_row_shape, pixels_per_chip_col_shape)

        # Index into chip_priorities, and add on the pixel priority
        chip_priorities[chip_row, chip_col] += pixel_priorities[index]


    return chip_priorities

def ChipsToCover(chip_layout):
    """
    Generates the chips to cover
    :param chip_layout: A numpy array of the proper chip shape
    :return: A linear list of (x,y) coordinates of nonzero chips
    """
    chips_to_cover = []
    for row_index, row in enumerate(chip_layout):
        for col_index, col in enumerate(row):
            if col > 0:
                chips_to_cover.append([row_index, col_index])
    return chips_to_cover

def PixelListToBitmap(pixel_list, pixel_priorities, num_pixels_row, num_pixels_col):
    bitmap_array = np.zeros((num_pixels_row,num_pixels_col))
    # Place ones int he positions where pixel_list indicates an activated pixel
    for index, pixel in enumerate(pixel_list):
        bitmap_array[pixel[0],pixel[1]] = pixel_priorities[index]
    return bitmap_array

def PlotBitmap(bitmap_array,num_pixels_row, num_pixels_col):
    plt.imshow(bitmap_array, interpolation="none", cmap="gray")
#    plt.xticks(np.arange(0,num_pixels_col))
#    plt.yticks(np.arange(0,num_pixels_row))
#    plt.show(block=False)

def BitmapToPixelList(filename):
    bitmap_image = Image.open(filename)
    bitmap_image = bitmap_image.convert('L') #TESTME
    pixel_numpy_array = np.array(bitmap_image)
    image_shape = np.shape(pixel_numpy_array)
    num_pixels_per_row = image_shape[0]
    num_pixels_per_col = image_shape[1]
    print "read %s as %d x %d pixels, rows x cols" % (filename, num_pixels_per_row, num_pixels_per_col)
# black is high priority, white is low
    pixel_priorities = 255 - pixel_numpy_array
    pixel_coordinate_list = []
    pixel_priorities_list = []
    rows,cols = np.nonzero(pixel_priorities)
    for i in range(0,len(rows)):
        pixel_coordinate_list.append((rows[i],cols[i]))
  	pixel_priorities_list.append( pixel_priorities[rows[i],cols[i]] )
#    print "pixel coordinate list = ", pixel_coordinate_list
#    print "pixel priority list = ", pixel_priorities_list
#    print "pixel priorities array = \n", pixel_priorities
# plot raw image 
    PlotBitmap(pixel_numpy_array, num_pixels_per_row, num_pixels_per_col)
    return pixel_coordinate_list, pixel_priorities_list, image_shape 


def VisualizeAllChipLayouts(pixel_list, pixel_priorities, pixel_shape, chip_mask_shape):

    pixel_rows = pixel_shape[0]
    pixel_cols = pixel_shape[1]

    # Convert to a numpy array

    bitmap = PixelListToBitmap(pixel_list, pixel_priorities, pixel_shape[0], pixel_shape[1])

    PlotBitmap(bitmap, pixel_rows, pixel_cols)

    left_shift, right_shift, down_shift, up_shift = FindMaxShifts(pixel_list, pixel_rows, pixel_cols)

    pixel_configuration_list = FindAllPixelConfigurations(pixel_list, left_shift, right_shift, down_shift,
                                                                      up_shift)

    chip_rows = chip_mask_shape[0]
    chip_cols = chip_mask_shape[1]


    pixels_per_chip_row = pixel_shape[0]/chip_mask_shape[0]
    pixels_per_chip_col = pixel_shape[1]/chip_mask_shape[1]
    for pix_list in pixel_configuration_list:
        print("The current pixel list is", pix_list)
        chips = FindActivatedChips(pix_list, pixels_per_chip_row, pixels_per_chip_col, chip_rows, chip_cols)
        print("The corresponding chips are", chips)
        VisualizeChipsAndPixels(pix_list, pixel_priorities, chips, pixel_shape[0], pixel_shape[1], chip_mask_shape[0], chip_mask_shape[1])
        #VisualizeActivatedChips(chips, chip_rows, chip_cols)


def VisualizeActivatedChips(chip_list, number_of_rows_of_chips, number_of_cols_of_chips):
    fig = plt.figure(1)
    plt.title("Subfootprints Placed")
    ax = fig.gca()
    ax.set_yticks(np.arange(0, number_of_rows_of_chips + 1, 1))
    ax.set_xticks(np.arange(0, number_of_cols_of_chips + 1, 1))
    for row_index, row in enumerate(chip_list):
        for col_index, col in enumerate(row):
            if col == 1:

                # Specify the x, y position of the chip
                # Remember, 0,0 is at the top left. M,0 is at the bottom left.
                x_position = col_index
                y_position = number_of_rows_of_chips - row_index - 1
                print("The x,y position is", x_position , " , ", y_position)
                ax.add_patch(Rectangle((x_position, y_position), 1, 1, facecolor="grey"))

    plt.grid()
#    plt.show(block=False)

def VisualizeChipsAndPixels(pixel_list, pixel_priorities, chip_list, number_of_rows_of_pixels, number_of_cols_of_pixels, number_of_rows_of_chips, number_of_cols_of_chips):
    fig = plt.figure(1)
    plt.title("Subfootprints Placed")
    ax = fig.gca()
    ax.set_yticks(np.arange(0, number_of_rows_of_chips + 1, 1))
    ax.set_xticks(np.arange(0, number_of_cols_of_chips + 1, 1))
    for row_index, row in enumerate(chip_list):
        for col_index, col in enumerate(row):
            if col == 1:

                # Specify the x, y position of the chip
                # Remember, 0,0 is at the top left. M,0 is at the bottom left.
                x_position = col_index
                y_position = number_of_rows_of_chips - row_index - 1
                print("The x,y position is", x_position , " , ", y_position)
                ax.add_patch(Rectangle((x_position, y_position), 1, 1, facecolor="grey"))

    plt.grid()

    fig2 = plt.figure(2)
    bitmap_array = PixelListToBitmap(pixel_list, pixel_priorities, number_of_rows_of_pixels, number_of_cols_of_pixels)

    plt.imshow(bitmap_array, interpolation="none", cmap="gray")
    plt.xticks(np.arange(0,number_of_cols_of_pixels))
    plt.yticks(np.arange(0,number_of_rows_of_pixels))
    plt.grid()
#    plt.show(block=False)

def SaveToBitmap(bitmap_array, filename_string):
    # Stupid PIL does not accept float data, so cast as uint8
    result = Image.fromarray(bitmap_array.astype(np.uint8))
    if ".bmp" not in filename_string:
        filename_string += ".bmp"
    result.save(filename_string)

