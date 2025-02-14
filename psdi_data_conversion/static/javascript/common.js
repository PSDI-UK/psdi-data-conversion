/*
  convert_common.js
  Created 2025/02/14 by Bryan Gillis

  Common functions used throughout the project
*/

export function getServiceMode() {
  // Gets whether or not the app is operating in "Service mode" - this is the mode used for the public web app.
  // Returns bool - True indicates service mode, False indicates local mode
  return sessionStorage.getItem("service_mode");
}