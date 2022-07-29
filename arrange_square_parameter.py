# Reevu Adakroy
# July 23, 2022
# arrange fibers in square orientation, approximated to match given parameters

from PIL import Image, ImageDraw
import math
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
# find W, H: W = circle count across, H = circle count down
def solve_WH(fiber_count, height, width):
    # Requirement 1:
        # W*H >= fiber_count
        # fiber_count > W*(H-1)
        # fiber_count > (W-1)*H
    # generate all integer pairs that satisfy Requirement 1:
    lst = WH_pairings(fiber_count)
    # Requirement 2: Minimize absolute difference between W/H and width/height
    ratio = max(width/height,height/width)
    min_difference = ratio
    min_index = 0
    for index in range(len(lst)):
        if abs(lst[index][1]/lst[index][0]-ratio)<min_difference:
            min_difference = abs(lst[index][1]/lst[index][0]-ratio)
            min_index = index

    if height>width: return lst[min_index][1],lst[min_index][0]
    return lst[min_index]

# return list of positive integer pairs a,b such that given k:
    # a*b >= k
    # k > a*(b-1)
    # k > (a-1)*b
    # a <= b
def WH_pairings(k):
    # try all a such that: 1 <= a < 0.5 + sqrt(k+0.25)
    return [(a,math.ceil(k/a)) for a in range(1, int(0.5+(k+0.25)**0.5)+1)]

def place_circles(H, W, n):
    return [(n/2+j*n, n/2+i*n) for i in range(H) for j in range(W)]

def main():
    fiber_count, volume_fraction, height, width = validate_input()
    H,W = solve_WH(fiber_count, height, width)
    radius = math.sqrt(height*width*volume_fraction/math.pi/fiber_count)
    n = min(width/W, height/H)

    # store placement as a list of circle center coordinates
    circle_lst = place_circles(H,W,n)
    if len(circle_lst)>fiber_count:
        print('It is recommended that your fiber count be adjusted to match ' + str(len(circle_lst)) + '.')
    else:
        print('Good choice of fiber count.')
    circle_lst = circle_lst[:int(fiber_count)]

    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))
    image.show()
    image.save('square.png')

if __name__ == '__main__': main()