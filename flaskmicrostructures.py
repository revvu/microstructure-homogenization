from flask import Flask, render_template, url_for, jsonify, request
import arrange_square, arrange_hexagonal, arrange_random, prepare_fvdam
import time

app = Flask(__name__)

@app.route("/")
def index():
    return render_template('index.html')


@app.route('/_generate_square')
def make_square():
    vol_frac = request.args.get('volumeFraction',0,type=float)
    time_saved = str(int(time.time()))
    filename = "generated_square"+time_saved+".png"
    with open('time_saved.txt', 'w') as file:
        print(time_saved, file=file)

    height, width = arrange_square.execute(vol_frac, time_saved)
    return jsonify(nothing=filename, pixel_height=height, pixel_width=width)

@app.route('/_generate_hexagonal')
def make_hexagonal():
    vol_frac = request.args.get('volumeFraction',0,type=float)
    time_saved = str(int(time.time()))
    filename = "generated_hexagonal"+time_saved+".png"

    with open('time_saved.txt', 'w') as file:
        print(time_saved, file=file)

    height, width = arrange_hexagonal.execute(vol_frac, time_saved)
    return jsonify(nothing=filename, pixel_height=height, pixel_width=width)

@app.route('/_generate_random')
def make_random():
    vol_frac = request.args.get('volumeFraction',0,type=float)
    fiber_count = request.args.get('fiberCount',0,type=int)
    height = request.args.get('height',0,type=int)
    width = request.args.get('width',0,type=int)


    time_saved = str(int(time.time()))
    filename = "generated_random"+time_saved+".png"
    with open('time_saved.txt', 'w') as file:
        print(time_saved, file=file)
    arrange_random.execute(vol_frac, fiber_count, height, width, time_saved)
    return jsonify(nothing=filename, pixel_height=height, pixel_width=width)

@app.route('/_prepare_fvdam')
def prep_fvdam():
    with open('time_saved.txt','r') as file:
        lines = file.readlines()
    time_saved = lines[0].strip()
    filename, approx_vol_frac = prepare_fvdam.main(time_saved)
    return jsonify(nothing=filename, approxVolumeFraction = str(approx_vol_frac))

# @app.route('/_run_fvdam')
# def execute_fvdam():
#     with open('time_saved.txt','r') as file:
#         lines = file.readlines()
#     time_saved = lines[0].strip()
#     LOP = request.args.get('loadingOption',type=int)
#     filename = run_fvdam.main(time_saved,LOP)
#     return jsonify(nothing=filename)


if __name__=='__main__':
    app.run(debug=True)
