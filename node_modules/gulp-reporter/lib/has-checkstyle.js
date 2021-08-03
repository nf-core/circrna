'use strict';
const path = require('path');
const fs = require('fs');
const jenkinsHome = process.env.JENKINS_HOME;
let hasCheckstyle;

try {
	hasCheckstyle = jenkinsHome && fs.statSync(path.join(jenkinsHome, 'plugins/checkstyle')).isDirectory();
} catch (ex) {
	//
}
module.exports = hasCheckstyle || false;
