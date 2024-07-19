/*
  feedback.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the feedback.htm gui work.
*/

$(document).ready(function() {
    const font = sessionStorage.getItem("font"),
          size = sessionStorage.getItem("size"),
          weight = sessionStorage.getItem("weight"),
          letter = sessionStorage.getItem("letter"),
          line = sessionStorage.getItem("line"),
          colour = sessionStorage.getItem("colour"),
          back = sessionStorage.getItem("back");

    if (font != null) {
        $(".normalText, .middle").css({
            fontFamily: font,
            fontSize: size,
            fontWeight: weight,
            letterSpacing: letter,
            lineHeight: line
        });

        $(".normalText").css({color: colour});
        $("h1, h2").css({letterSpacing: letter});
        $("h1").css({color: colour});
        $("form, #upper").css({background: back});
    }
});

