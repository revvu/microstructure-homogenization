# Reevu Adakroy
# July 25, 2022
# arrange fibers in square orientation with chosen parameters and inputted volume fraction

from PIL import Image, ImageDraw
import math
import sys; args = sys.argv[1:]

# make input checking more robust
def validate_input():
    if len(args) != 1:
        print('Invalid input.')
        exit()
    return float(args[0])

# find W, H: W = circle count across, H = circle count down
def solve_WH(fiber_count, height, width):
    lst = WH_pairings(fiber_count)
    ratio = max(width/height,height/width)
    min_difference = ratio
    min_index = 0
    for index in range(len(lst)):
        if abs(lst[index][1]/lst[index][0]-ratio)<min_difference:
            min_difference = abs(lst[index][1]/lst[index][0]-ratio)
            min_index = index

    if height>width: return lst[min_index][1],lst[min_index][0]
    return lst[min_index]

def WH_pairings(k):
    return [(a,math.ceil(k/a)) for a in range(1, int(0.5+(k+0.25)**0.5)+1)]

def place_circles(H, W, n):
    return [(n/2+j*n, n/2+i*n) for i in range(H) for j in range(W)]

def main():
    volume_fraction = validate_input()
    fiber_count = 12
    height = 200
    width = 150
    H,W = solve_WH(fiber_count, height, width)
    radius = math.sqrt(height*width*volume_fraction/math.pi/fiber_count)
    n = min(width/W, height/H)

    # store placement as a list of circle center coordinates
    circle_lst = place_circles(H,W,n)
    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))
    image.show()
    image.save('square.png')

def execute(vol_frac, time_saved):
    volume_fraction = vol_frac
    fiber_count = 12
    height = 150
    width = 200
    H,W = solve_WH(fiber_count, height, width)
    radius = math.sqrt(height*width*volume_fraction/math.pi/fiber_count)
    n = min(width/W, height/H)

    # store placement as a list of circle center coordinates
    circle_lst = place_circles(H,W,n)
    image = Image.new('RGBA', (width, height))

    draw = ImageDraw.Draw(image)
    draw.rectangle((0,0,width,height),fill=(54, 79, 107))
    for circle in circle_lst:
        draw.ellipse((circle[0] - radius, circle[1] - radius, circle[0] + radius, circle[1] + radius), fill=(63, 193, 201))

    filename = 'generated_square'+str(time_saved)+'.png'
    image.save('static/images/'+filename)

    with open('static/Generated Microstructures/RUC/generated_microstructure_'+time_saved+'.txt', 'w') as file:
        print(height,file=file)
        print(width,file=file)
        print(radius,file=file)
        for circle in circle_lst:
            print(circle[0],circle[1],file=file)

    return height, width

if __name__ == '__main__': main()
