/*
  feedback.js
  Version 1.0, 26th June 2024

  This is the JavaScript which makes the feedback.htm gui work.
*/

import { loadAccessibilitySettings } from './accessibility.js';

const r = document.querySelector(':root');

$(document).ready(function() {
    loadAccessibilitySettings();
});

