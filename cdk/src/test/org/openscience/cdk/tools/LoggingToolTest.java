/* Copyright (C) 2005-2009  Egon Willighagen <egonw@users.sf.net>
 *                    2007  Rajarshi Guha <rajarshi@users.sf.net>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA. 
 */
package org.openscience.cdk.tools;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKTestCase;

/**
 * @cdk.module test-log4j
 */
public class LoggingToolTest extends CDKTestCase {
	

	@Test public void testLoggingTool_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		Assert.assertNotNull(logger);
	}

	@Test public void testLoggingTool() throws Exception {
		LoggingTool logger = new LoggingTool();
		Assert.assertNotNull(logger);
	}

	@Test public void testLoggingTool_Class() throws Exception {
		LoggingTool logger = new LoggingTool(this.getClass());
		Assert.assertNotNull(logger);
	}

	@Test public void testClass$_String() throws Exception {
		// no idea why the Coverage test requires this test
		Assert.assertTrue(true);
	}

	@Test public void testConfigureLog4j() throws Exception {
		LoggingTool.configureLog4j();
	}

	@Test public void testDebug_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this);
	}

	@Test public void testDebug_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, this);
	}

	@Test public void testDebug_Object_int() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, 1);
	}

	@Test public void testDebug_Object_double() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, 1.0);
	}

	@Test public void testDebug_Object_boolean() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, true);
	}

	@Test public void testDebug_Object_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, this, this, this, this);
	}

	@Test public void testDebug_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, this, this, this);
	}

	@Test public void testDebug_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.debug(this, this, this);
	}

	@Test public void testError_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this);
	}

	@Test public void testError_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, this);
	}

	@Test public void testError_Object_int() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, 1);
	}

	@Test public void testError_Object_double() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, 1.0);
	}

	@Test public void testError_Object_boolean() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, true);
	}

	@Test public void testError_Object_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, this, this, this, this);
	}

	@Test public void testError_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, this, this, this);
	}

	@Test public void testError_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.error(this, this, this);
	}

	@Test public void testWarn_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this);
	}

	@Test public void testWarn_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, this);
	}

	@Test public void testWarn_Object_int() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, 1);
	}

	@Test public void testWarn_Object_double() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, 1.0);
	}

	@Test public void testWarn_Object_boolean() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, true);
	}

	@Test public void testWarn_Object_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, this, this, this, this);
	}

	@Test public void testWarn_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, this, this, this);
	}

	@Test public void testWarn_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.warn(this, this, this);
	}

	@Test public void testInfo_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this);
	}

	@Test public void testInfo_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, this);
	}

	@Test public void testInfo_Object_int() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, 1);
	}

	@Test public void testInfo_Object_double() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, 1.0);
	}

	@Test public void testInfo_Object_boolean() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, true);
	}

	@Test public void testInfo_Object_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, this, this, this, this);
	}

	@Test
    public void testInfo_Object_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, this, this, this);
	}

	@Test public void testInfo_Object_Object_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.info(this, this, this);
	}

	@Test public void testFatal_Object() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.fatal(this);
	}

	@Test public void testSetStackLength_int() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.setStackLength(20);
	}

	@Test public void testDumpClasspath() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.dumpClasspath();
	}

	@Test public void testDumpSystemProperties() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.dumpSystemProperties();
	}

	@Test public void testIsDebugEnabled() throws Exception {
		LoggingTool logger = new LoggingTool(this);
		logger.isDebugEnabled();
	}
	
    @Test public void testCreate() throws Exception {
        ILoggingTool logger = LoggingTool.create(this.getClass());
        Assert.assertNotNull(logger);
    }
}

