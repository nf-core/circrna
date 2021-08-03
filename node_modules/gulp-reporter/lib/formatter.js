'use strict';
const chalk = require('chalk');
const path = require('path');
const emphasize = require('emphasize');
const cliTruncate = require('cli-truncate');
const stringWidth = require('string-width');
const ansiEscapes = require('ansi-escapes');
const isCI = require('ci-info').isCI;
const icons = require('./icons');
let fixIcon;

if (!isCI && process.stdout.isTTY && require('is-windows')()) {
	if (process.env.VSCODE_PID || process.env.TERM_PROGRAM === 'Hyper') {
		fixIcon = str => {
			if (str === icons.error) {
				str = str[0] + ' ';
			} else if (str !== icons.info) {
				str += ' ';
			}
			return str;
		};
	} else if (process.env.ConEmuPID) {
		fixIcon = str => {
			if (str === icons.info) {
				str = '‼';
			}
			return str;
		};
	} else {
		fixIcon = str => {
			str = str[0];
			return str;
		};
	}
} else {
	fixIcon = str => {
		if (str === icons.error) {
			str = str[0];
		}
		return str;
	};
}

const msgColor = {
	success: chalk.green,
	info: chalk.blue,
	warn: chalk.yellow,
	error: chalk.red,
};

let padNum;
if (String.prototype.padStart) {
	padNum = (number, targetLength) => {
		return String(number).padStart(targetLength, '0');
	};
} else {
	padNum = (number, targetLength) => {
		number = String(number);
		return ('0').repeat(Math.max(targetLength - number.length, 0)) + number;
	};
}

function justify (left, right, width) {
	if (left) {
		width -= stringWidth(left);
	} else {
		left = '';
	}
	if (right) {
		width -= stringWidth(right);
	} else {
		right = '';
	}
	return left + ' '.repeat(Math.max(width, left && right ? 1 : 0)) + right;
}

function timeFormat (timestamp) {
	timestamp = new Date(timestamp * 1000);
	return [
		timestamp.getFullYear(),
		padNum(timestamp.getMonth() + 1, 2),
		padNum(timestamp.getDate(), 2),
	].join('-') + ' ' + [
		padNum(timestamp.getHours(), 2),
		padNum(timestamp.getMinutes(), 2),
		padNum(timestamp.getSeconds(), 2),
	].join(':');
}

function highlight (error) {
	let source = error.source.replace(/\t/g, '    ');
	if (chalk.enabled) {
		const extname = path.extname(error.fileName).slice(1);
		if (extname) {
			try {
				return emphasize.highlight(extname.slice(1), source).value;
			} catch (ex) {
				//
			}
		}
		source = emphasize.highlightAuto(source).value;
	}
	return source;
}

function addLineNumbers (code, line, fmtOpts) {
	function add (s) {
		let prefix = ' '.repeat(fmtOpts.columnNumberWidth + 1);
		prefix += chalk.yellow(padNum(line, fmtOpts.lineNumberWidth) + ' | ');
		line++;
		return prefix + s;
	}
	return code.replace(/^.*$/gm, add);
}

function formatError (error, fmtOpts) {
	const message = [];
	let severity = error.severity;
	if (error.demote || !severity || !icons[severity]) {
		severity = 'success';
	}

	const mainMessage = chalk.gray(`${
		padNum(error.lineNumber || 1, fmtOpts.lineNumberWidth)
	}:${
		padNum(error.columnNumber || 1, fmtOpts.columnNumberWidth)
	}`) + ' ' + msgColor[severity](icons[severity] + ' ' + error.message);
	let subMsg = [];

	if (error.plugin) {
		subMsg.push(justify(error.plugin, null, fmtOpts.maxPluginNameWidth));
	}

	if (error.rule) {
		subMsg.push(justify(error.rule, null, fmtOpts.maxRuleNameWidth));
	}

	if (subMsg.length === 1) {
		subMsg[0] = justify(subMsg[0], null, fmtOpts.maxPluginNameWidth + fmtOpts.maxRuleNameWidth + 1);
	}

	const link = error.docShort || error.doc;

	if (link) {
		subMsg.push(chalk.underline(link));
	}
	if (subMsg.length) {
		subMsg = `(${subMsg.join(' ').trim()})`;
		if (stringWidth(mainMessage) + stringWidth(subMsg) >= fmtOpts.termColumns) {
			message.push(mainMessage, justify(
				null,
				subMsg,
				fmtOpts.termColumns
			));
		} else {
			message.push(
				justify(
					mainMessage,
					subMsg,
					fmtOpts.termColumns
				)
			);
		}
	} else {
		message.push(mainMessage);
	}

	message[0] = message[0].replace(/\S\uFE0F/g, fixIcon);

	if (fmtOpts.source && error.source) {
		let source = fmtOpts.sourceCache[error.source];
		if (!source) {
			source = highlight(error);
			source = addLineNumbers(source, error.lineNumber || 1, fmtOpts);
			source = cliTruncate(source, fmtOpts.termColumns);
			fmtOpts.sourceCache[error.source] = source;
		}
		message.push(source);
	}

	if (error.blame) {
		let blame = error.blame._string;
		if (!blame) {
			const author = `(${
				justify(
					error.blame.author.name,
					null,
					fmtOpts.maxNameWidth
				)
			} ${
				justify(
					`<${error.blame.author.mail}>`,
					null,
					fmtOpts.maxEmailWidth + 2
				)
			} ${
				timeFormat(error.blame.author.time)
			})`;

			blame = chalk.gray(justify(
				cliTruncate(
					error.blame.hash,
					Math.max(fmtOpts.termColumns - stringWidth(author) - 1, 9)
				),
				author,
				fmtOpts.termColumns
			));
			error.blame._string = blame;
		}
		message.push(blame);
	}

	return message.join('\n');
}

function formatter (file, options) {
	const sourceCache = {};
	let maxLineNumber = 0;
	let maxColumnNumber = 0;
	let maxNameWidth = 0;
	let maxEmailWidth = 0;
	let maxPluginNameWidth = 0;
	let maxRuleNameWidth = 0;
	let maxMessageWidth = 0;
	let errors = file.report.errors;
	errors.forEach(error => {
		if (error.lineNumber > maxLineNumber) {
			maxLineNumber = error.lineNumber;
		}
		if (error.columnNumber > maxColumnNumber) {
			maxColumnNumber = error.columnNumber;
		}
		if (error.blame && !/^0+$/.test(error.blame.hash)) {
			const nameWidth = stringWidth(error.blame.author.name);
			if (nameWidth > maxNameWidth) {
				maxNameWidth = nameWidth;
			}
			const emailWidth = stringWidth(error.blame.author.mail);
			if (emailWidth > maxEmailWidth) {
				maxEmailWidth = emailWidth;
			}
		}
		const pluginNameWidth = error.plugin && stringWidth(error.plugin);
		if (pluginNameWidth > maxPluginNameWidth) {
			maxPluginNameWidth = pluginNameWidth;
		}
		const ruleNameWidth = error.rule && stringWidth(error.rule);
		if (ruleNameWidth > maxRuleNameWidth) {
			maxRuleNameWidth = ruleNameWidth;
		}
		const messageWidth = stringWidth(error.message);
		if (messageWidth > maxMessageWidth) {
			maxMessageWidth = messageWidth;
		}
	});
	const lineNumberWidth = Math.max(String(maxLineNumber).length, 2);
	const columnNumberWidth = Math.max(String(maxColumnNumber).length, 2);
	const fmtOpts = {
		sourceCache,
		maxNameWidth,
		maxEmailWidth,
		maxPluginNameWidth,
		maxRuleNameWidth,
		lineNumberWidth,
		columnNumberWidth,
		termColumns: Math.min(
			options._termColumns - 6,
			Math.max(
				31 + lineNumberWidth + columnNumberWidth + maxMessageWidth + maxPluginNameWidth + maxRuleNameWidth,
				options.blame ? 66 + maxNameWidth + maxEmailWidth : 0
			)
		),
		source: options.blame,
	};

	errors = errors.map(error => (
		formatError(error, fmtOpts).replace(/^/gm, '    ')
	));

	let head = chalk.cyan.underline(path.relative(file.cwd, file.path).replace(/\\/g, '/'));
	if (process.stdout.isTTY && !isCI) {
		// make relative paths Cmd+click'able in iTerm
		head = ansiEscapes.iTerm.setCwd(file.cwd) + head;
	}
	errors.unshift(head);
	return errors.join('\n');
}

module.exports = formatter;
