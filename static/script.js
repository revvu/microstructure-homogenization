function squareScript() {
  $(".fvdam").css("display","none");
  $(".microstructure-info").css("display","block");
  $(".square-info").css("display","block");
  $(".hexagonal-info").css("display","none");
  $(".random-info").css("display","none");
  $(".microstructure-img").css("transform","none");
  $("#square-img").css("transform","scale(1.25)");

  $("#square-img-result").css("display", "block");
  $("#hexagonal-img-result").css("display", "none");
  $("#random-img-result").css("display", "none");
}

function hexagonalScript() {
  $(".fvdam").css("display","none");
  $(".microstructure-info").css("display","block");
  $(".square-info").css("display","none");
  $(".hexagonal-info").css("display","block");
  $(".random-info").css("display","none");
  $(".microstructure-img").css("transform","none");
  $("#hexagonal-img").css("transform","scale(1.25)");

  $("#square-img-result").css("display", "none");
  $("#hexagonal-img-result").css("display", "block");
  $("#random-img-result").css("display", "none");
}

function randomScript() {
  $(".fvdam").css("display","none");
  $(".microstructure-info").css("display","block");
  $(".square-info").css("display","none");
  $(".hexagonal-info").css("display","none");
  $(".random-info").css("display","block");
  $(".microstructure-img").css("transform","none");
  $("#random-img").css("transform","scale(1.25)");

  $("#square-img-result").css("display", "none");
  $("#hexagonal-img-result").css("display", "none");
  $("#random-img-result").css("display", "block");
}

function update_square(){
  $.getJSON($SCRIPT_ROOT + '/_generate_square', {
    volumeFraction: $('input[name="volumeFraction"]').val()
  }, function(data) {
    $(".img-result").attr("src","static/images/"+data.nothing);
    $(".img-result").height(2*data.pixel_height);
    $(".img-result").width(2*data.pixel_width);
  });
  $(".fvdam").css("display","block");
  return false;
}

function update_hexagonal(){
  $.getJSON($SCRIPT_ROOT + '/_generate_hexagonal', {
    volumeFraction: $('input[name="volumeFraction"]').val()
  }, function(data) {
    $(".img-result").attr("src","static/images/"+data.nothing);
    $(".img-result").height(2*data.pixel_height);
    $(".img-result").width(2*data.pixel_width);
  });
  $(".fvdam").css("display","block");
  return false;
}

function update_random(){
  $(".generate-loading").show();
  $(".fvdam").css("display","block");
  $(".prep-results").css("display","none");
  $.getJSON($SCRIPT_ROOT + '/_generate_random', {
    volumeFraction: $('input[name="volumeFraction"]').val(),
    fiberCount: $('input[name="fiberCount"]').val(),
    height: $('input[name="height"]').val(),
    width: $('input[name="width"]').val(),
  }, function(data) {
    $(".img-result").attr("src","static/images/"+data.nothing);
    $(".img-result").height(2*data.pixel_height);
    $(".img-result").width(2*data.pixel_width);
    $(".generate-loading").hide();
  });

  return false;
}

function prep_fvdam(){
  $.getJSON($SCRIPT_ROOT + '/_prepare_fvdam', {

  }, function(data) {
    $(".prep-results").attr("href",data.nothing);
  });
  $(".prep-results").css("display","inline-block");
  return false;
}

function execute_fvdam(){
  $(".delay-run").text("Running...");
  $(".final-results").css("display","none");
  $(".loading").show();
  $.getJSON($SCRIPT_ROOT + '/_run_fvdam', {
     loadingOption: $('select[name="loading-option"]').val()
  }, function(data) {
    $(".final-results").attr("href",data.nothing);
    $(".delay-run").text("Finished.");
    $(".final-results").css("display","inline-block");
    $(".loading").hide();
  });

  return false;
}
