/*
  feedback.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the feedback.htm gui work.
*/

import { connectModeToggleButton } from './accessibility.js';

$(document).ready(function() {
  connectModeToggleButton();
});

// Remove the loading cover when everything is loaded
$(window).on('load', function() {
    $("#cover").hide();
});

