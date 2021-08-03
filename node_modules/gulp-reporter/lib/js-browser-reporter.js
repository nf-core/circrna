'use strict';
/* eslint-env browser */
/* eslint no-console: "off" */
/* eslint no-var: "off" */

/**
 * 代码错误汇报函数，在浏览器中运行，用于将收集到的js的错误信息在浏览器控制台中弹出
 * 注意！！此函数将被toString()后发送到浏览器，并非在node下运行！！
 * @param  {LintError[]}  errors  错误信息
 * @return {void}
 */
function jsBrowserReporter (errors) {
	var isGecko = !!window.netscape;

	var throwError = isGecko ? function (message, stack, err) {
		stack = stack.map((uri) => {
			return '@' + uri;
		});
		if (!err.severity || err.severity === 'error') {
			var error = new (window.SyntaxError || Error)(message);
			error.columnNumber = err.columnNumber;
			error.fileName = err.fileName;
			error.lineNumber = err.lineNumber;
			error.stack = stack.join('\n');
			setTimeout(() => {
				throw error;
			}, 0);
		} else {
			consoleError(message, stack, err);
		}
	} : function (message, stack, err) {
		stack = stack.map((uri) => {
			return '    at ' + uri;
		});
		consoleError(message, stack, err);
	};

	function consoleError (message, stack, err) {
		stack.unshift(err.name + ': ' + message);
		stack = stack.join('\n');
		var level = err.demote ? 'log' : (err.severity in console ? err.severity : 'error');
		console[level](stack);
	}

	// 将文件路径与模块路径拼接为完整的url
	errors.forEach((err) => {
		var message = err.message + ' (' + [
			err.plugin,
			err.rule,
		].filter(Boolean).join(' ') + ')';
		var stack = [
			[
				'file://' + err.fileName.replace(/^(\w:)\\/, '/$1/').replace(/\\/g, '/'),
				err.lineNumber,
				err.columnNumber,
			].filter(Boolean).join(':'),
		];
		if (err.doc) {
			stack.push(err.doc);
		}

		throwError(message, stack, err);
	});
}

module.exports = jsBrowserReporter.toString().replace(/^(function)\s*\w+/, '$1');
