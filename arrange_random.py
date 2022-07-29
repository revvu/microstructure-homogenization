# Reevu Adakroy
# July 25, 2022
# arrange fibers in random orientation, starting with square approximation then perturbing

from PIL import Image, ImageDraw
import arrange_square_parameter
import math
import random
import sys; args = sys.argv[1:]

# make input checking more robust
def validate_input():
    if len(args) != 4:
        print('Invalid input.')
        exit()
    lst = [float(arg) for arg in args]
    for i in range(len(lst)):
        if i!=1: lst[i] = int(lst[i])
    return lst

def set_speeds(circle_lst):
    multiple = 1
    return [(2*multiple*random.random()-multiple, 2*multiple*random.random()-multiple) for i in range(len(circle_lst))]

def perturb(circle_lst, speeds, height, width, radius):
    buffered_radius = 1.1*radius
    # check pairwise circles and bounce apart when needed
    for i in range(len(circle_lst)):
        for j in range(i+1, len(circle_lst)):
            if (circle_lst[i][0]-circle_lst[j][0])**2+(circle_lst[i][1]-circle_lst[j][1])**2<=4*buffered_radius**2:
                if (circle_lst[j][0]-circle_lst[i][0])/(speeds[j][0]-speeds[i][0])<0 or (circle_lst[j][1]-circle_lst[i][1])/(speeds[j][1]-speeds[i][1])<0:
                    temp = speeds[i]
                    speeds[i] = speeds[j]
                    speeds[j] = temp

    for i in range(len(circle_lst)):
        circle_lst[i] = (circle_lst[i][0]+speeds[i][0],circle_lst[i][1] + speeds[i][1])

        # check for collision with left wall
        if circle_lst[i][0] - buffered_radius < 0:
            speeds[i] = (abs(speeds[i][0]), speeds[i][1])
            circle_lst[i] = (circle_lst[i][0] + 2 * speeds[i][0], circle_lst[i][1])

        # check for collision with upper wall
        if circle_lst[i][1] - buffered_radius < 0:
            speeds[i] = (speeds[i][0], abs(speeds[i][1]))
            circle_lst[i] = (circle_lst[i][0], circle_lst[i][1] + 2 * speeds[i][1])

        # check for collision with right wall
        if circle_lst[i][0] + buffered_radius > width:
            speeds[i] = (-abs(speeds[i][0]), speeds[i][1])
            circle_lst[i] = (circle_lst[i][0] + 2 * speeds[i][0], circle_lst[i][1])

        # check for collision with bottom wall
        if circle_lst[i][1] + buffered_radius > height:
            speeds[i] = (speeds[i][0], -abs(speeds[i][1]))
            circle_lst[i] = (circle_lst[i][0], circle_lst[i][1] + 2 * speeds[i][1])

def check_overlap(circle_lst, radius):
    buffered_radius = 1.1 * radius
    for i in range(len(circle_lst)):
        for j in range(i+1, len(circle_lst)):
            if (circle_lst[i][0]-circle_lst[j][0])**2+(circle_lst[i][1]-circle_lst[j][1])**2<=4*buffered_radius**2: return True
    return False

def main():
    fiber_count, volume_fraction, height, width = validate_input()
    H,W = arrange_square_parameter.solve_WH(fiber_count, height, width)
    radius = math.sqrt(height*width*volume_fraction/math.pi/fiber_count)
    n = min(width/W, height/H)

    # store placement as a list of circle center coordinates
    circle_lst = arrange_square_parameter.place_circles(H,W,n)[:int(fiber_count)]
    speeds = set_speeds(circle_lst)

    for i in range(1000): perturb(circle_lst,speeds,height, width, radius)

    # while there is overlap, keep perturbing
    while check_overlap(circle_lst, radius): perturb(circle_lst,speeds,height, width, radius)

    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))

    image.show()
    image.save('random.png')


def execute(volume_fraction, fiber_count, height, width, time_saved):
    H,W = arrange_square_parameter.solve_WH(fiber_count, height, width)
    radius = math.sqrt(height*width*volume_fraction/math.pi/fiber_count)
    n = min(width/W, height/H)

    # store placement as a list of circle center coordinates
    circle_lst = arrange_square_parameter.place_circles(H,W,n)[:int(fiber_count)]
    speeds = set_speeds(circle_lst)

    for i in range(1000): perturb(circle_lst,speeds,height, width, radius)

    # while there is overlap, keep perturbing
    while check_overlap(circle_lst, radius): perturb(circle_lst,speeds,height, width, radius)

    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))

    image.save('static/images/generated_random'+time_saved+'.png')

    with open('static/Generated Microstructures/RUC/generated_microstructure_'+time_saved+'.txt', 'w') as file:
        print(height,file=file)
        print(width,file=file)
        print(radius,file=file)
        for circle in circle_lst:
            print(circle[0],circle[1],file=file)


if __name__ == '__main__': main()