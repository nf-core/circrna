'use strict';
const LintError = require('./lint-error');

/**
 *
 *
 * @class TSLintError
 * @extends {LintError}
 */
class TSLintError extends LintError {
	constructor (error) {
		super({
			// 文件名
			// fileName: error.fileName,

			// TSLint无警告，错误等级全部算错误
			severity: 'error',

			// 行号
			lineNumber: error.startPosition.lineAndCharacter.line + 1,

			// 列号
			columnNumber: error.startPosition.lineAndCharacter.character + 1,

			// 错误信息
			message: error.failure,

			// 触发错误的规则
			rule: error.ruleName,

			// 源代码上下文
			source: error.rawLines,

			// 报错插件
			plugin: 'TSLint',

			// 文档
			doc: error.ruleName && `https://palantir.github.io/tslint/rules/${error.ruleName}/`,
		}, error);
	}
}

module.exports = TSLintError;
